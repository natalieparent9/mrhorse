

# JAGS model

mvmr_horse_model_jags = function() {
  for (i in 1:N){
    by[i] ~ dnorm(mu[i], 1 / (sy[i] * sy[i]))
    mu[i] = inprod(bx0[i, 1:K], theta) + alpha[i]
    bx[i,1:K] ~ dmnorm(bx0[i,1:K], Tx[1:K, ((i-1)*K+1):(i*K)])

    kappa[i] = (rho[i]^2 / (1 + K*rho[i]^2))
    bx0[i,1:K] ~ dmnorm(mx + sx0 * rho[i] * alpha[i] / (phi[i] * tau), A - kappa[i] * B)
    r[i] ~ dbeta(10, 10);T(, 1)
    rho[i] = 2*r[i] -1
    alpha[i] ~ dnorm(0, 1 / (tau * tau * phi[i] * phi[i]))
    phi[i] = a[i] / sqrt(b[i])
    a[i] ~ dnorm(0, 1);T(0, )
    b[i] ~ dgamma(0.5, 0.5)
  }

  c ~ dnorm(0, 1);T(0, )
  d ~ dgamma(0.5, 0.5)
  tau = c / sqrt(d)

  mx ~ dmnorm(rep(0, K), R[,])

  for (k in 1:K){
    vx0[k] ~ dnorm(0, 1);T(0, )
    sx0[k] = sqrt(vx0[k])
    theta[k] ~ dunif(-10, 10)
    for (j in 1:K){
      A[j, k] = ifelse(j==k, 1/vx0[j], 0)
      B[j, k] = 1 / (sx0[j] * sx0[k])
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
#' # Check estimates
#' # MREx$MR_Estimate
#'
#' # Fit model with Stan
#' # MREx = mvmr_horse(data_mv_ex, n.warmup = 1000, n.iter = 2000, stan = TRUE)
mvmr_horse = function(D, n.chains = 3, variable.names = "theta", n.iter = 10000, n.burnin = 10000, stan = FALSE, n.cores = parallelly::availableCores()){

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
  if (!all(c('betaY','betaYse') %in% colnames(D))) stop("Error: betaY and/or betaYse are missing")

  # Ensure at least the results for theta parameter are saved
  if("theta" %in% variable.names){
    variable.names = variable.names
  } else{
    variable.names = c("theta", variable.names)
  }

  # Initialize precision matrix
  Tx = matrix(0, nrow = K, ncol = p * K)
  for (j in 1:p) {
    Tx[, ((j-1)*K+1):(j*K)] = diag(Sx[j, ]^2)
  }

  data_list = data = list(by = D$betaY, bx = Bx, sy = D$betaYse, Tx = Tx, N = p, K = K, R = diag(K))

  cat("Fitting model with ", K, " exposures", sep='')

  if (stan == TRUE) {  # Run using Stan
    fit = rstan::sampling(stanmodels$mvmr_horse,
                        data = data_list,
                        pars = variable.names,
                        chains = n.chains,
                        iter = n.burnin + n.iter,
                        warmup = n.burnin,
                        cores = n.cores,
                        n.thin = 1                  # remove thinning that occurs by default
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
                                       n.cluster = n.cores)
      mr.coda = coda::as.mcmc(fit)

  } else {
    # stop("Error: invalid input for argument stan, must be TRUE or FALSE")
  }


  mr_estimate = data.frame("Parameter" = sprintf("theta[%i]", 1:K),
                           "Estimate" = unname(summary(mr.coda)$statistics[sprintf("theta[%i]", 1:K), 1]),
                           "SD" = unname(summary(mr.coda)$statistics[sprintf("theta[%i]", 1:K), 2]),
                           "2.5% quantile" = unname(summary(mr.coda)$quantiles[sprintf("theta[%i]", 1:K), 1]),
                           "97.5% quantile" = unname(summary(mr.coda)$quantiles[sprintf("theta[%i]", 1:K), 5]),
                           "Rhat" = unname(coda::gelman.diag(mr.coda)$psrf[sprintf("theta[%i]", 1:K), 1]))
  names(mr_estimate) = c("Parameter", "Estimate", "SD", "2.5% quantile", "97.5% quantile", "Rhat")
  mr_estimate[,-1] = round(mr_estimate[,-1], 3)

  return(list('fit' = fit, "MR_Estimate" = mr_estimate, "MR_Coda" = mr.coda))
}


