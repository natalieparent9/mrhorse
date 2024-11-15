

# JAGS model

mvmr_horse_model_jags = function() {
  for (i in 1:N){                                             # Number of variants (rows)
    for (k in 1:K) {
      mu[i,k] = bx0[i,k]                                      # Expected mean of exposures (bxs)

      # Create an array of N K+1xK+1 matrices, each variant has a K+1xK+1 matrix
      for (j in 1:K) {                                        # Fills KxK diagonal with variance of bxs and off diagonal with 0
        cov[i, k, j] = ifelse(k == j, sx[i,k] * sx[i,k], 0)
      }
      cov[i, k, K+1] =  omega * sx[i,k] * sy[i]               # Covariance between bxs and by - last column of matrix
      cov[i, K+1, k] =  cov[i, k, K+1]                        # Covariance between bxs and by - last row of matrix
    }
    cov[i, K+1, K+1] = sy[i] * sy[i]                          # Inserts var of by of by into last diagonal element
    mu[i, K+1] = inprod(bx0[i, 1:K], theta) + alpha[i]        # Expected mean of by

    obs[i,] ~ dmnorm.vcov(mu[i,], cov[i,,])                   # Jointly model bxs and by


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

  tau = ifelse(fixed_tau == -1, c / sqrt(d), fixed_tau)    # Fixed or estimated depending on user specification                                                            # Controls the global level of shrinkage for alphas

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


#' Multivariable MR-Horse method
#'
#' @param D Dataset containing betaY, betaYse and at least betaX1, betaX1se, betaX2 and betaX2se, or MRMVInput object
#' @param n.chains Number of chains
#' @param variable.names Parameters to save estimates for, in addition to theta
#' @param n.iter Number of iterations (not including warmup)
#' @param n.burnin Number of warmup iterations
#' @param stan fit the model using stan, default is to use JAGS
#' @param n.cores Number of cores to use in parallel if running multiple chains, defaults to parallelly::availableCores()
#' @param return_fit return the fitted JAGS or Stan model object, default is FALSE
#' @param fixed_tau a fixed value for tau, otherwise tau will be estimated
#' @param omega correlation parameter for bx and by, default is 0 which assumes independence
#'
#' @return Output from the mr_horse() function is a list that contains:
#' $MR_Estimate: a data frame with the causal effect estimate (which is the posterior mean), standard deviation (i.e., the posterior standard deviation), upper and lower bounds of the 95% credible interval, and the R-hat value                                                                                                               the posterior standard deviation), upper and lower bounds of the 95% credible interval, and the R-hat value
#  $MR_Coda: full MCMC samples for all parameters in `variable.names`
#  $Fit: the fitted model object (if reuturn_fit was set to TRUE)
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
mvmr_horse = function(D, n.chains = 3, variable.names = "theta", n.iter = 10000, n.burnin = 10000, stan = FALSE, n.cores = parallelly::availableCores(), return_fit = FALSE, fixed_tau = -1, omega=0){

  # Validate input
  if (is.data.frame(D)) {                        # accepts data.frame, tibble, data.table
    N = dim(D)[1]                                # number of variants (rows)
    K = sum(grepl("^betaX[0-9]+$", names(D)))    # number of exposures

    Bx = D[, sprintf("betaX%i", 1:K)]            # dataframe of associations with exposures
    Sx = D[, sprintf("betaX%ise", 1:K)]          # dataframe of standard errors of associations with exposures


  } else if (inherits(D, "MRMVInput")) {         # If D is an MRMVInput object
    N = length(D@betaY)                          # number of variants (rows)
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

  data_list = list(obs = as.matrix(D[, c(sprintf("betaX%i", 1:K), 'betaY')], ncol=2),
                   sx = Sx,
                   sy = D$betaYse,
                   N = N,
                   K = K,
                   R = diag(K),
                   fixed_tau=fixed_tau,
                   omega=omega)

  cat("Fitting model with ", K, " exposures and ", N, " variants\n", sep='')

  if (stan == TRUE) {  # Run using Stan
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
    fit = R2jags::jags.parallel(data = data_list,
                                     model.file = mvmr_horse_model_jags,
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
                           "Rhat" = unname(coda::gelman.diag(mr.coda[,sprintf("theta[%i]", 1:K)])$psrf[1:K]))
  names(mr_estimate) = c("Parameter", "Estimate", "SD", "2.5% quantile", "97.5% quantile", "Rhat")
  mr_estimate[,-1] = round(mr_estimate[,-1], 3)

  if (return_fit == TRUE) return(list("MR_Estimate" = mr_estimate, "MR_Coda" = mr.coda, 'Fit' = fit))
  return(list("MR_Estimate" = mr_estimate, "MR_Coda" = mr.coda))
}


