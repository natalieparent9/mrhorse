library(R2jags)
library(rstan)
library(coda)
library(parallelly)


mr_horse_model_jags = function() {
  for (i in 1:N){                                # Where N is the number of variants (rows)
    by[i] ~ dnorm(mu[i], 1/(sy[i] * sy[i]))      # Likelihood for by
    mu[i] = theta * bx0[i] + alpha[i]            # Mean for bx0 distribution
    bx[i] ~ dnorm(bx0[i], 1 / (sx[i] * sx[i]))   # Likelihood for bx

    bx0[i] ~ dnorm(mx0 + (sqrt(vx0)/(tau * phi[i])) * rho[i] * alpha[i], 1 / ((1 - rho[i]^2) * vx0)) # Effect of the variant on the exposure
    r[i] ~ dbeta(10, 10);T(, 1)                  # Correlation between alpha and bx0
    rho[i] = 2*r[i] -1                           # Converts the truncated beta distribution parameter r[i] to a correlation parameter rho[i].

    alpha[i] ~ dnorm(0, 1 / (tau * tau * phi[i] * phi[i]))  # Pleiotropic effects
    phi[i] = a[i] / sqrt(b[i])                   # Scaling factor for each variant based on parameters a and b
    a[i] ~ dnorm(0, 1);T(0, )                    # Parameter for phi[i]
    b[i] ~ dgamma(0.5, 0.5)                      # Parameter for phi[i]
  }

  c ~ dnorm(0, 1);T(0, )       # Parameter for tau
  d ~ dgamma(0.5, 0.5)         # Parameter for tau
  tau = c / sqrt(d)            # Controls the global level of shrinkage for alphas

  vx0 ~ dnorm(0, 1);T(0, )     # Variance for bx0 distribution
  mx0 ~ dnorm(0, 1)            # Mean for bx0 distribution

  theta ~ dunif(-10, 10)       # Causal effect of X on Y
}


# Stan model can be found in /inst/stan


#' Univariable MR-Horse method
#'
#' Causal effect estimation is performed in a Bayesian framework using a horseshoe shrinkage prior to account for both correlated and uncorrelated pleiotropic effects
#'
#' @param D Data frame containing betaY, betaYse, betaX and betaXse, or MRInput object
#' @param n.chains Number of chains
#' @param variable.names Parameters to save estimates for
#' @param n.iter Number of iterations (not including warmup)
#' @param n.burnin Number of warmup iterations
#' @param stan Fit the model using stan, default is to use JAGS
#' @param n.cores Number of cores to use in parallel if running multiple chains, defaults to parallelly::availableCores()
#'
#' @return List
#' @export
#'
#' @examples
#' library(coda)
#'
#' # Load example data
#' data(data_ex)
#'
#' # Fit model with JAGS
#' MREx = mr_horse(data_ex)
#'
#' # Check estimates
#' MREx$MR_Estimate
mr_horse = function(D, n.chains = 3, variable.names = "theta", n.iter = 10000, n.burnin = 10000, stan = FALSE, n.cores = parallelly::availableCores()){

  # Validate input
  if (is.data.frame(D)) {  # accepts data.frame, tibble, data.table
  } else if (inherits(D, "MRInput")) {
    # If D is an MRInput object, convert it to a data frame
    D = data.frame(
      betaX = D@betaX,
      betaY = D@betaY,
      betaXse = D@betaXse,
      betaYse = D@betaYse
    )
  } else {
    stop("Error: D must be a data frame or an MRInput object.")
  }

  if (all(c('betaX', 'betaY', 'betaXse', 'betaYse') %in% colnames(D))) {
  } else {
    stop("Error: D must contain columns: betaX, betaY, betaXse, betaYse")
  }

  # Ensure at least the results for theta parameter are saved
  if("theta" %in% variable.names){
    variable.names = variable.names
  } else{
    variable.names = c("theta", variable.names)
  }

  # Fit model
  if (stan == TRUE) {
    fit = rstan::sampling(stanmodels$mr_horse,
                   pars = variable.names,
                   data = list("N"=nrow(D), by=D$betaY, bx = D$betaX, sy = D$betaYse, sx = D$betaXse),
                   iter = n.iter + n.burnin,
                   warmup = n.burnin,
                   chains = n.chains,
                   # control = list(adapt_delta = 0.9, max_treedepth = 12),
                   cores = n.cores)
    mr.coda = rstan::As.mcmc.list(fit)

  } else if (stan == FALSE) {
      fit = R2jags::jags.parallel(data = list(by = D$betaY, bx = D$betaX, sy = D$betaYse, sx = D$betaXse, N = length(D$betaY)),
                  parameters.to.save = variable.names,
                  n.chains = n.chains,
                  n.iter = n.burnin + n.iter,
                  n.burnin = n.burnin,
                  model.file = mr_horse_model_jags,
                  n.cluster = n.cores)
      mr.coda = coda::as.mcmc(fit)
  } else {
      stop("Error: invalid input for argument stan, must be TRUE or FALSE")
  }

  mr_estimate = data.frame("Estimate" = summary(mr.coda[, "theta"])$statistics[[1]],
                           "SD" = summary(mr.coda[, "theta"])$statistics[[2]],
                           "2.5% quantile"  = summary(mr.coda[, "theta"])$quantiles[[1]],
                           "97.5% quantile" = summary(mr.coda[, "theta"])$quantiles[[5]],
                           "Rhat" = coda::gelman.diag(mr.coda)$psrf[[1]])
  mr_estimate = round(mr_estimate, 3)
  names(mr_estimate) = c("Estimate", "SD", "2.5% quantile", "97.5% quantile", "Rhat")
  return(list("MR_Estimate" = mr_estimate, "MR_Coda" = mr.coda))
}
