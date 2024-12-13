
# JAGS model
# Note: can use ifelse() in JAGS model but not regular if() statements

mr_horse_model_jags = function() {
    for (i in 1:J){                                  # Where J is the number of variants (rows)

      mu[i, 1] = bx0[i]                              # Expected mean of bx distribution
      mu[i, 2] = theta * bx0[i] + alpha[i]           # Expected mean of by distribution

      obs[i, ] ~ dmnorm.vcov(mu[i,], V[,,i])         # Jointly model bx and by likelihood, .vcov specifies we have provided var covar matrix instead of precision

      bx0[i] ~ dnorm(mx0 + (sqrt(vx0)/(tau * phi[i])) * rho[i] * alpha[i], 1 / ((1 - rho[i]^2) * vx0)) # Effect of the variant on the exposure
      r[i] ~ dbeta(10, 10);T(0, 1)                   # Will be scaled to correlation between alpha and bx0
      rho[i] = 2*r[i] -1                             # Converts the truncated beta distribution parameter r[i] to a correlation parameter rho[i].

      alpha[i] ~ dnorm(0, 1 / (tau * tau * phi[i] * phi[i]))  # Pleiotropic effects
      phi[i] = a[i] / sqrt(b[i])                     # Scaling factor for each variant based on parameters a and b
      a[i] ~ dnorm(0, 1);T(0, )                      # Parameter for phi[i]
      b[i] ~ dgamma(0.5, 0.5)                        # Parameter for phi[i]
    }

    c ~ dnorm(0, 1);T(0, )                           # Parameter for tau
    d ~ dgamma(0.5, 0.5)                             # Parameter for tau

    tau = ifelse(fixed_tau == -1, c / sqrt(d), fixed_tau)    # Fixed or estimated depending on user specification

    vx0 ~ dnorm(0, 1);T(0, )                         # Variance for bx0 distribution
    mx0 ~ dnorm(0, 1)                                # Mean for bx0 distribution

    theta ~ dunif(-10, 10)                           # Causal effect of X on Y
}


# Stan model can be found in /inst/stan


#' Univariable MR-Horse method
#'
#' Causal effect estimation is performed in a Bayesian framework using a horseshoe shrinkage prior to account for both correlated and uncorrelated pleiotropic effects
#'
#' @param D Data frame containing betaY, betaYse, betaX and betaXse, or MRInput object
#' @param n.chains Number of chains
#' @param variable.names Parameters to save estimates for, in addition to theta
#' @param n.iter Number of iterations (not including warmup)
#' @param n.burnin Number of warmup iterations
#' @param stan Fit the model using stan, default is to use JAGS
#' @param n.cores Number of cores to use in parallel if running multiple chains, defaults to parallelly::availableCores()
#' @param return_fit Return the fitted JAGS or Stan model object, default is FALSE
#' @param fixed_tau Optional fixed value for tau, otherwise tau will be estimated
#' @param omega Optional correlation parameter for betaX and betaY, default is 0 which assumes independence
#'
#' @return Output from the mr_horse() function is a list that contains:
#' $MR_Estimate: a data frame with the causal effect estimate (which is the posterior mean), standard deviation (i.e., the posterior standard deviation), upper and lower bounds of the 95% credible interval, and the R-hat value                                                                                                               the posterior standard deviation), upper and lower bounds of the 95% credible interval, and the R-hat value
#  $MR_Coda: full MCMC samples for all parameters in `variable.names`
#  $Fit: the fitted model object (if reuturn_fit was set to TRUE)
#' @export
#'
#' @examples
#' # Fit model with JAGS
#' # data(data_ex)
#' # MREx = mr_horse(data_ex, n.burnin = 1000, n.iter = 2000)
#'
#' # Check estimates
#' # MREx$MR_Estimate
#'
#' # Fit model with Stan
#' # MREx = mr_horse(data_ex, n.warmup = 1000, n.iter = 2000, stan = TRUE)
#'
#' # Save samples for tau parameter (in addition to theta)
#' # MREx = mr_horse(data_ex, n.warmup = 1000, n.iter = 2000, variable.names = 'tau')
#' # summary(MREx$MR_Coda)
#'
#' # View diagnostic plots (trace and density)
#' # plot(MREx$MR_Coda[,c('theta','tau')])
#'
mr_horse = function(D, n.chains = 3, variable.names = "theta", n.iter = 10000, n.burnin = 10000,
                    stan = FALSE, n.cores = parallelly::availableCores(), return_fit = FALSE, fixed_tau = -1, omega = 0){

  J = nrow(D) # Number of genetic instruments

  # Validate data input, accepts data frame, data table, tibble or MRInput object (does not require conversion)
  if (!(is.data.frame(D) | inherits(D, "MRInput"))) stop("Error: D must be a data frame or an MRInput object")

  # Convert MRInput to data frame
  if (inherits(D, "MRInput")) D = data.frame(betaX = D@betaX, betaY = D@betaY, betaXse = D@betaXse, betaYse = D@betaYse)

  # Check variables and arguments
  vars = c('betaX', 'betaY', 'betaXse', 'betaYse')
  if (!(all(vars %in% colnames(D)))) stop("Error: D must contain columns: betaX, betaY, betaXse, betaYse")
  if (fixed_tau != -1 & fixed_tau < 0) stop('Error: fixed value for tau must be >0')
  if (!(omega >= -1 && omega <= 1)) stop('Error: omega must be between -1 and 1')

  # Ensure at least the results for the theta parameter are saved
  variable.names = unique(c("theta", variable.names))

  # Create covariance matrix V
  R = matrix(c(1, omega, omega, 1), nrow=2,ncol=2)      # R matrix with omega on off diagonals, default is identity
  V = array(0, c(2,2,J))
  for (i in 1:J) {
    S = diag(D[i, c('betaXse', 'betaYse')], nrow = 2, ncol=2)
    V[,,i] = S %*% R %*% S
  }
  l = list()
  for (i in 1:J) l[[i]]= V[,,i]                          # Stan requires Nx2x2 format rather than 2x2xN

  data_list = list(
    J = J,                                               # Number of genetic instruments
    obs = as.matrix(D[,c('betaX','betaY')],ncol=2),      # Observed bx and by
    V = if (stan == TRUE) l else V,                      # Covariance matrix of betaX and betaY
    fixed_tau = fixed_tau                                # Fixed tau value - default -1
  )

  # Fit model
  if (stan == TRUE) {
    fit = rstan::sampling(stanmodels$mr_horse,
                   pars = variable.names,
                   data = data_list,
                   iter = n.iter + n.burnin,
                   warmup = n.burnin,
                   chains = n.chains,
                   open_progress = FALSE, # prevent pop up during testing, does not suppress progress being printed in console
                   # control = list(adapt_delta = 0.9, max_treedepth = 12),
                   cores = n.cores)
    mr.coda = rstan::As.mcmc.list(fit)

  } else if (stan == FALSE) {
      fit = R2jags::jags.parallel(data = data_list,
                  parameters.to.save = variable.names,
                  n.chains = n.chains,
                  n.iter = n.burnin + n.iter,
                  n.burnin = n.burnin,
                  n.thin = 1, # retain all samples instead of thinning to 10%
                  model.file = mr_horse_model_jags,
                  n.cluster = n.cores)
      mr.coda = coda::as.mcmc(fit)
  } else {
      stop("Error: invalid input for argument stan, must be TRUE or FALSE")
  }

  mr_estimate = data.frame("Estimate"       = summary(mr.coda[, "theta"])$statistics[[1]],
                           "SD"             = summary(mr.coda[, "theta"])$statistics[[2]],
                           "2.5% quantile"  = summary(mr.coda[, "theta"])$quantiles[[1]],
                           "97.5% quantile" = summary(mr.coda[, "theta"])$quantiles[[5]],
                           "Rhat"           = coda::gelman.diag(mr.coda[, "theta"])$psrf[[1]])
  mr_estimate = round(mr_estimate, 3)
  names(mr_estimate) = c("Estimate", "SD", "2.5% quantile", "97.5% quantile", "Rhat")

  if (return_fit == TRUE) return(list("MR_Estimate" = mr_estimate, "MR_Coda" = mr.coda, 'Fit' = fit))
  return(list("MR_Estimate" = mr_estimate, "MR_Coda" = mr.coda))
}
