
# JAGS model
# The function is set up so that different versions of the model can be used depending on options selected
# Code corresponding to the macro chosen will be subsituted into the model file, since JAGS models do not allow conditionals within

mr_horse_model_jags = function(fixed_tau, omega) {

  model = "function() {
    for (i in 1:N){                                # Where N is the number of variants (rows)

      mu[i, 1] = bx0[i]                              # Expected mean of bx distribution
      mu[i, 2] = theta * bx0[i] + alpha[i]           # Expected mean of by distribution
      vary[i] = sy[i] * sy[i]                        # Variance of by
      varx[i] = sx[i] * sx[i]                        # Variance of bx
      cov[i] = $OMEGA * sx[i] * sy[i]                # Covariance of bx and by
      precision[i, 1, 1] = vary[i] / (varx[i] * vary[i] - cov[i]^2)    # Precision matrix (inverse of covariance matrix)
      precision[i, 2, 2] = varx[i] / (varx[i] * vary[i] - cov[i]^2)
      precision[i, 1, 2] = -cov[i] / (varx[i] * vary[i] - cov[i]^2)
      precision[i, 2, 1] = precision[i, 1, 2]
      obs[i, ] ~ dmnorm(mu[i,], precision[i,,])      # Jointly model bx and by likelihood

      bx0[i] ~ dnorm(mx0 + (sqrt(vx0)/(tau * phi[i])) * rho[i] * alpha[i], 1 / ((1 - rho[i]^2) * vx0)) # Effect of the variant on the exposure
      r[i] ~ dbeta(10, 10);T(, 1)                    # Correlation between alpha and bx0
      rho[i] = 2*r[i] -1                             # Converts the truncated beta distribution parameter r[i] to a correlation parameter rho[i].

      alpha[i] ~ dnorm(0, 1 / (tau * tau * phi[i] * phi[i]))  # Pleiotropic effects
      phi[i] = a[i] / sqrt(b[i])                     # Scaling factor for each variant based on parameters a and b
      a[i] ~ dnorm(0, 1);T(0, )                      # Parameter for phi[i]
      b[i] ~ dgamma(0.5, 0.5)                        # Parameter for phi[i]
    }

    c ~ dnorm(0, 1);T(0, )       # Parameter for tau
    d ~ dgamma(0.5, 0.5)         # Parameter for tau
    $TAU                         # Fixed or estimated depending on user specification

    vx0 ~ dnorm(0, 1);T(0, )     # Variance for bx0 distribution
    mx0 ~ dnorm(0, 1)            # Mean for bx0 distribution

    theta ~ dunif(-10, 10)       # Causal effect of X on Y
  }"

  # Sub in options
  tau_model = if (fixed_tau == -1) { # uses argument passed to mr_horse function
    'tau = c / sqrt(d)'
  } else {
    'tau = fixed_tau' # the model will find the fixed_tau value as its included in the data list
  }

  model = gsub("$TAU", tau_model, model, fixed=TRUE)
  model = gsub("$OMEGA", omega, model, fixed=TRUE)   # uses argument passed to mr_horse function
  model
}


# Stan model can be found in /inst/stan


#' Univariable MR-Horse method
#'
#' Causal effect estimation is performed in a Bayesian framework using a horseshoe shrinkage prior to account for both correlated and uncorrelated pleiotropic effects
#'
#' @param D data frame containing betaY, betaYse, betaX and betaXse, or MRInput object
#' @param n.chains number of chains
#' @param variable.names parameters to save estimates for, in addition to theta
#' @param n.iter number of iterations (not including warmup)
#' @param n.burnin number of warmup iterations
#' @param stan fit the model using stan, default is to use JAGS
#' @param n.cores number of cores to use in parallel if running multiple chains, defaults to parallelly::availableCores()
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
#' #
#'
mr_horse = function(D, n.chains = 3, variable.names = "theta", n.iter = 10000, n.burnin = 10000,
                    stan = FALSE, n.cores = parallelly::availableCores(), return_fit = FALSE, fixed_tau = -1, omega = 0){

  # Validate data input, accepts data frame, data table, tibble or MRInput object (does not require conversion)
  if (!(is.data.frame(D) | inherits(D, "MRInput"))) stop("Error: D must be a data frame or an MRInput object")

  vars = c('betaX', 'betaY', 'betaXse', 'betaYse')
  if (!(all(vars %in% colnames(D)) | all(vars %in% slotNames(D)))) stop("Error: D must contain columns: betaX, betaY, betaXse, betaYse")
  if (fixed_tau != -1 & fixed_tau < 0) stop('Error: fixed value for tau must be >0')
  if (!(omega >= -1 && omega <= 1)) stop('Error: omega must be between -1 and 1')

  # Ensure at least the results for theta parameter are saved
  variable.names = unique(c("theta", variable.names))

  data_list = list(N=length(D$betaY), by=D$betaY, bx = D$betaX, sy = D$betaYse, sx = D$betaXse, fixed_tau = fixed_tau,
                   obs = as.matrix(D[,c('betaX','betaY')], ncol=2))

  # Fit model
  if (stan == TRUE) {
    fit = rstan::sampling(
      stanmodels$mr_horse,
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
      if (fixed_tau == -1) data_list$fixed_tau = NULL # avoid warning of unused variable
      fit = R2jags::jags.parallel(data = data_list,
                  parameters.to.save = variable.names,
                  n.chains = n.chains,
                  n.iter = n.burnin + n.iter,
                  n.burnin = n.burnin,
                  model.file = eval(parse(text=mr_horse_model_jags(fixed_tau, omega))),
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
