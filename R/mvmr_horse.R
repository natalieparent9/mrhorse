#' library(R2jags)
#' library(rstan)
#' library(coda)
#' library(parallelly)
#'
#'
#' mvmr_horse_model_jags = function() {
#'   for (i in 1:N){
#'     by[i] ~ dnorm(mu[i], 1 / (sy[i] * sy[i]))
#'     mu[i] = inprod(bx0[i, 1:K], theta) + alpha[i]
#'     bx[i,1:K] ~ dmnorm(bx0[i,1:K], Tx[1:K, ((i-1)*K+1):(i*K)])
#'
#'     kappa[i] = (rho[i]^2 / (1 + K*rho[i]^2))
#'     bx0[i,1:K] ~ dmnorm(mx + sx0 * rho[i] * alpha[i] / (phi[i] * tau), A - kappa[i] * B)
#'     r[i] ~ dbeta(10, 10);T(, 1)
#'     rho[i] = 2*r[i] -1
#'     alpha[i] ~ dnorm(0, 1 / (tau * tau * phi[i] * phi[i]))
#'     phi[i] = a[i] / sqrt(b[i])
#'     a[i] ~ dnorm(0, 1);T(0, )
#'     b[i] ~ dgamma(0.5, 0.5)
#'   }
#'
#'   c ~ dnorm(0, 1);T(0, )
#'   d ~ dgamma(0.5, 0.5)
#'   tau = c / sqrt(d)
#'
#'   mx ~ dmnorm(rep(0, K), R[,])
#'
#'   for (k in 1:K){
#'     vx0[k] ~ dnorm(0, 1);T(0, )
#'     sx0[k] = sqrt(vx0[k])
#'     theta[k] ~ dunif(-10, 10)
#'     for (j in 1:K){
#'       A[j, k] = ifelse(j==k, 1/vx0[j], 0)
#'       B[j, k] = 1 / (sx0[j] * sx0[k])
#'     }
#'   }
#' }
#'
#'
#' #' Title
#' #'
#' #' @param D Dataset containing betaY, betaYse and at least betaX1, betaX1se, betaX2 and betaX2se
#' #' @param n.chains Number of chains
#' #' @param variable.names Parameters to save estimates for
#' #' @param n.iter Number of iterations (not including warmup)
#' #' @param n.burnin Number of warmup iterations
#' #' @param n.chains Number of cores to use in parallel if running multiple chains, defaults to parallelly::availableCores()
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' mvmr_horse = function(D, n.chains = 3, variable.names = "theta", n.iter = 10000, n.burnin = 10000, n.cores = parallelly::availableCores()){
#'
#'
#'   if("theta" %in% variable.names){
#'     variable.names = variable.names
#'   } else{
#'     variable.names = c("theta", variable.names)
#'   }
#'
#'   p = dim(D)[1]
#'   K = sum(sapply(1:dim(D)[2], function(j){substr(names(D)[j], 1, 5)=="betaX"}))/2
#'
#'   Bx = D[, sprintf("betaX%i", 1:K)]
#'   Sx = D[, sprintf("betaX%ise", 1:K)]
#'   Tx = matrix(nrow = K, ncol = p*K)
#'   for (j in 1:p){
#'     Tx[, ((j-1)*K+1):(j*K)] = diag(1 / Sx[j, ]^2)
#'   }
#'
#'   jags_fit = R2jags::jags.parallel(data = list(by = D$betaY, bx = Bx, sy = D$betaYse, Tx = Tx, N = p, K = K, R = diag(K)),
#'                           parameters.to.save = variable.names,
#'                           n.chains = no_ini,
#'                           n.iter = n.burnin + n.iter,
#'                           n.burnin = n.burnin,
#'                           model.file = mvmr_horse_model_jags,
#'                           n.cluster = n.chains)
#'   mr.coda = coda::as.mcmc(jags_fit)
#'   mr_estimate = data.frame("Parameter" = sprintf("theta[%i]", 1:K),
#'                            "Estimate" = round(unname(summary(mr.coda)$statistics[sprintf("theta[%i]", 1:K), 1]), 3),
#'                            "SD" = round(unname(summary(mr.coda)$statistics[sprintf("theta[%i]", 1:K), 2]), 3),
#'                            "2.5% quantile" = round(unname(summary(mr.coda)$quantiles[sprintf("theta[%i]", 1:K), 1]), 3),
#'                            "97.5% quantile" = round(unname(summary(mr.coda)$quantiles[sprintf("theta[%i]", 1:K), 5]), 3),
#'                            "Rhat" = round(unname(gelman.diag(mr.coda)$psrf[sprintf("theta[%i]", 1:K), 1]), 3))
#'   names(mr_estimate) = c("Parameter", "Estimate", "SD", "2.5% quantile", "97.5% quantile", "Rhat")
#'   return(list("MR_Estimate" = mr_estimate, "MR_Coda" = mr.coda))
#' }
