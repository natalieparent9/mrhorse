

# JAGS model

mvmr_horse_model_jags = function(f=c('standard', 'fixed_tau')) {
  f = match.arg(f)

  model = "function() {
    for (i in 1:N){                                                   # number of variants (rows)
      by[i] ~ dnorm(mu[i], 1 / (sy[i] * sy[i]))                       # Association between the exposures and outcome
      mu[i] = inprod(bx0[i, 1:K], theta) + alpha[i]                   # Mean effect on Y for each variant
      bx[i,1:K] ~ dmnorm(bx0[i,1:K], Tx[1:K, ((i-1)*K+1):(i*K)])      # Association between the variants and exposures

      kappa[i] = (rho[i]^2 / (1 + K*rho[i]^2))                        # Used to adjust bx0
      bx0[i,1:K] ~ dmnorm(mx + sx0 * rho[i] * alpha[i] / (phi[i] * tau), A - kappa[i] * B) # Latent effect of the variants on the exposures
      r[i] ~ dbeta(10, 10);T(, 1)                                     # Correlation between alpha and bx0
      rho[i] = 2*r[i] -1                                              # Converts the truncated beta distribution parameter r[i] to a correlation parameter rho[i]
      alpha[i] ~ dnorm(0, 1 / (tau * tau * phi[i] * phi[i]))          # Pleiotropic effects on the outcome
      phi[i] = a[i] / sqrt(b[i])                                      # Scaling factor for each variant based on parameters a and b
      a[i] ~ dnorm(0, 1);T(0, )                                       # Parameter for phi[i]
      b[i] ~ dgamma(0.5, 0.5)                                         # Parameter for phi[i]
    }

    c ~ dnorm(0, 1);T(0, )                                            # Parameter for tau
    d ~ dgamma(0.5, 0.5)                                              # Parameter for tau
    $TAU                                                              # Controls the global level of shrinkage for alphas

    mx ~ dmnorm(rep(0, K), R[,])

    for (k in 1:K){
      vx0[k] ~ dnorm(0, 1);T(0, )              # Variance for bx0 distribution
      sx0[k] = sqrt(vx0[k])                    # Standard error for bx0 distribution
      theta[k] ~ dunif(-10, 10)                # Causal effect of X[k] on Y
      for (j in 1:K){
      A[j, k] = ifelse(j==k, 1/vx0[j], 0)    # Precision matrix
      B[j, k] = 1 / (sx0[j] * sx0[k])        # Covariance matrix
      }
    }
  }
  "
  macros = list(list("$TAU",
               switch(f,
                      standard='tau = c / sqrt(d)',
                      fixed_tau='tau = fixed_tau')))

  for (m in seq(macros)) {
    model = gsub(macros[[m]][1], macros[[m]][2], model, fixed=TRUE) # macros[[m]][1] is the name e.g. 'standard' and macros[[m]][2] is corresponding code
  }
  model
}


#' Multivariable MR-Horse method
#'
#' @param D Dataset containing betaY, betaYse and at least betaX1, betaX1se, betaX2 and betaX2se, or MRMVInput object
#' @param n.chains Number of chains
#' @param variable.names Parameters to save estimates for, in addition to theta
#' @param n.iter Number of iterations (not including warmup)
#' @param n.burnin Number of warmup iterations
#' @param stan fit the model using stan, default is to use JAGS
#' @param n.cores Number of cores to use in parallel if running multiple chains, defaults to parallelly::availableCores()
#' @param fixed_tau a fixed value for tau, otherwise tau will be estimated
#' @param return_fit return the fitted JAGS or Stan model object, default is FALSE
#'
#' @return Output from the mr_horse() function is a list that contains:
#' $MR_Estimate: a data frame with the causal effect estimate (which is the posterior mean), standard deviation (i.e., the posterior standard deviation), upper and lower bounds of the 95% credible interval, and the R-hat value                                                                                                               the posterior standard deviation), upper and lower bounds of the 95% credible interval, and the R-hat value
#  $MR_Coda: full MCMC samples for all parameters in `variable.names`
#' @export
#'
#' @examples
#' # Fit model with JAGS
#' # data(data_mv_ex)
#' # MREx = mvmr_horse(data_mv_ex, n.warmup = 1000, n.iter = 2000)
#'
#' # Check estimates for theta
#' # MREx$MR_Estimate
#'
#' # Fit model with Stan
#' # MREx = mvmr_horse(data_mv_ex, n.warmup = 1000, n.iter = 2000, stan = TRUE)
#'
#' # View diagnostic plots (trace and density)
#' # plot(MREx$MR_Coda[,c('theta[1]','theta[2]')])
#'
mvmr_horse = function(D, n.chains = 3, variable.names = "theta", n.iter = 10000, n.burnin = 10000, stan = FALSE, n.cores = parallelly::availableCores(), fixed_tau = -1, return_fit = FALSE){

  # Validate input
  if (is.data.frame(D)) {                        # accepts data.frame, tibble, data.table
    p = dim(D)[1]                                # number of variants (rows)
    K = sum(grepl("^betaX[0-9]+$", names(D)))    # number of exposures

    Bx = D[, sprintf("betaX%i", 1:K)]            # dataframe of associations with exposures
    Sx = D[, sprintf("betaX%ise", 1:K)]          # dataframe of standard errors of associations with exposures

  } else if (inherits(D, "MRMVInput")) {         # If D is an MRMVInput object
    p = length(D@betaY)                          # number of variants (rows)
    K = length(D@betaX[,1])                      # number of exposures

    Bx = as.data.frame(D@betaX[1:K,1])           # dataframe of associations with exposures
    colnames(Bx) = paste0("betaX", 1:K)
    Sx = as.data.frame(D@betaXse[1:K,1])         # dataframe of standard errors of associations with exposures
    colnames(Sx) = paste0("betaX", 1:K, 'se')

  } else {
    stop("Error: D must be a data frame or an MRMVInput object.")
  }

  if (K < 2) stop("Error: Dataset must contain at least two exposures")
  if (!all(c('betaY','betaYse') %in% colnames(D)) & !all(c('betaY','betaYse') %in% slotNames(D))) {
    stop("Error: betaY and/or betaYse are missing")
  }

  if (fixed_tau != -1 & fixed_tau < 0) stop('Error: fixed value for tau must be >0')

  # Ensure at least the results for theta parameter are saved
  variable.names = unique(c("theta", variable.names))

  # Initialize precision/covariance matrix
  Tx = matrix(0, nrow = K, ncol = p * K)

  data_list = list(by = D$betaY, bx = Bx, sy = D$betaYse, Tx = Tx, N = p, K = K, R = diag(K))

  # Handling of fixed or estimated tau parameter
  if (fixed_tau == -1) {  # estimate tau
    model = 'standard'    # for use in macros to control JAGS model
  } else {                # fix tau
    model = 'fixed_tau'
    data_list$fixed_tau = fixed_tau
  }

  cat("Fitting model with ", K, " exposures and ",p, " variants\n", sep='')

  if (stan == TRUE) {  # Run using Stan
    for (j in 1:p) {data_list$Tx[, ((j-1)*K+1):(j*K)] = diag(Sx[j, ]^2)} # fill covariance matrix
    fit = rstan::sampling(stanmodels$mvmr_horse,
                        data = data_list,
                        pars = variable.names,
                        chains = n.chains,
                        iter = n.burnin + n.iter,
                        warmup = n.burnin,
                        cores = n.cores,
                        open_progress = FALSE # prevent pop up during testing, does not suppress progress being printed in console
                        # control = ifelse(algorithm=='NUTS',list(adapt_delta = adapt_delta, max_treedepth = max_treedepth), list())
    )
    mr.coda = rstan::As.mcmc.list(fit)

  } else if (stan == FALSE) {  # Run using JAGS
    for (j in 1:p) {data_list$Tx[, ((j-1)*K+1):(j*K)] = diag(1 / Sx[j, ]^2)}   # fill precision matrix
    fit = R2jags::jags.parallel(data = data_list,
                                     model.file = eval(parse(text=mvmr_horse_model_jags(model))),
                                     parameters.to.save = variable.names,
                                     n.chains = n.chains,
                                     n.iter = n.burnin + n.iter,
                                     n.burnin = n.burnin,
                                     n.cluster = n.cores,
                                     n.thin = 1)                # remove thinning that occurs by default
    mr.coda = coda::as.mcmc(fit)

  } else {
      stop("Error: invalid input for argument stan, must be TRUE or FALSE")
  }


  mr_estimate = data.frame("Parameter" = sprintf("theta[%i]", 1:K),
                           "Estimate" = unname(summary(mr.coda)$statistics[sprintf("theta[%i]", 1:K), 1]),
                           "SD" = unname(summary(mr.coda)$statistics[sprintf("theta[%i]", 1:K), 2]),
                           "2.5% quantile" = unname(summary(mr.coda)$quantiles[sprintf("theta[%i]", 1:K), 1]),
                           "97.5% quantile" = unname(summary(mr.coda)$quantiles[sprintf("theta[%i]", 1:K), 5]),
                           "Rhat" = unname(coda::gelman.diag(mr.coda)$psrf[sprintf("theta[%i]", 1:K), 1]))
  names(mr_estimate) = c("Parameter", "Estimate", "SD", "2.5% quantile", "97.5% quantile", "Rhat")
  mr_estimate[,-1] = round(mr_estimate[,-1], 3)

  if (return_fit == TRUE) return(list("MR_Estimate" = mr_estimate, "MR_Coda" = mr.coda, 'Fit' = fit))
  return(list("MR_Estimate" = mr_estimate, "MR_Coda" = mr.coda))
}


