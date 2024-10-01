library(R2jags)
library(rstan)
library(coda)

mr_horse_model_jags = function() {
  for (i in 1:N){
    by[i] ~ dnorm(mu[i], 1/(sy[i] * sy[i]))
    mu[i] = theta * bx0[i] + alpha[i]
    bx[i] ~ dnorm(bx0[i], 1 / (sx[i] * sx[i]))

    bx0[i] ~ dnorm(mx0 + (sqrt(vx0)/(tau * phi[i])) * rho[i] * alpha[i], 1 / ((1 - rho[i]^2) * vx0))
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

  vx0 ~ dnorm(0, 1);T(0, )
  mx0 ~ dnorm(0, 1)

  theta ~ dunif(-10, 10)
}


mr_horse_model_stan = "
data {
  int<lower=0> N;        // Number of observations (SNPs)
  vector[N] by;          // Observed effect estimates on Y (beta_y)
  vector[N] bx;          // Observed effect estimates on X (beta_x)
  vector[N] sy;          // Standard errors for by
  vector[N] sx;          // Standard errors for bx
}

parameters {
  real theta;            // Causal effect
  vector[N] alpha;       // Pleiotropic effects
  vector[N] bx0;         // Latent effect estimates on X
  vector<lower=0>[N] a;  // Parameter for phi[i]
  vector<lower=0>[N] b;  // Parameter for phi[i]
  vector<lower=0, upper=1>[N] r;  // Truncated beta distribution parameter for correlation
  real<lower=0> c;       // Parameter for tau
  real<lower=0> d;       // Parameter for tau
  real<lower=0> vx0;     // Variance for bx0 distribution
  real mx0;              // Mean for bx0 distribution
}

transformed parameters {
  vector<lower=0>[N] phi = a ./ sqrt(b);        // Scaling factor for each SNP based on parameters a and b
  vector[N] rho = 2 * r - 1;           // Converts the truncated beta distribution parameter r[i] to a correlation parameter rho[i].
  vector[N] mu = theta * bx0 + alpha;  // Computes the mean effect on Y for each SNP
  real<lower=0> tau = c / sqrt(d);              // Precision parameter

}

model {
  // Priors
  theta ~ uniform(-10, 10);                    // Uniform prior for causal effect theta
  alpha ~ normal(0, tau * phi);                // Normal prior for alpha
  bx0 ~ normal(mx0 + (sqrt(vx0) / (tau * phi)) .* rho .* alpha, sqrt((1 - rho .* rho) * vx0)); // Normal prior for bx0
  a ~ normal(0, 1);                            // Normal prior for a, truncated to be positive
  b ~ gamma(0.5, 0.5);                         // Gamma prior for b
  r ~ beta(10, 10);                            // Beta prior for r, truncated to [0, 1] in params section
  c ~ normal(0, 1);                            // Normal prior for c, truncated to be positive
  d ~ gamma(0.5, 0.5);                         // Gamma prior for d
  vx0 ~ normal(0, 1);                          // Normal prior for vx0, truncated to be positive
  mx0 ~ normal(0, 1);                          // Normal prior for mx0

  // Likelihood
  by ~ normal(mu, sy);                         // Likelihood for by (observed beta_y)
  bx ~ normal(bx0, sx);                        // Likelihood for bx (observed beta_x)
}
"

#' Title
#'
#' @param D Data frame containing betaY, betaYse, betaX and betaXse, or MRInput object
#' @param no_ini Number of chains
#' @param variable.names Parameters to save estimates for
#' @param n.iter Number of iterations (not including warmup)
#' @param n.burnin Number of warmup iterations
#' @param stan Fit the model using stan, default is to use JAGS
#'
#' @return
#' @export
#'
#' @examples
mr_horse = function(D, no_ini = 3, variable.names = "theta", n.iter = 10000, n.burnin = 10000, stan = FALSE){
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

  # Ensure at least the results for theta are saved
  if("theta" %in% variable.names){
    variable.names = variable.names
  } else{
    variable.names = c("theta", variable.names)
  }







  if (stan == TRUE) {
    print('stan')
  } else if (stan == FALSE) {
    print('jags')
      # jags_fit = R2jags::jags(data = list(by = D$betaY, bx = D$betaX, sy = D$betaYse, sx = D$betaXse, N = length(D$betaY)),
      #             parameters.to.save = variable.names,
      #             n.chains = no_ini,
      #             n.iter = n.burnin + n.iter,
      #             n.burnin = n.burnin,
      #             model.file = mr_horse_model)
  } else {
      stop("Error: invalid input for argument stan, must be TRUE or FALSE")
  }



  # mr.coda = coda::as.mcmc(jags_fit)
  # mr_estimate = data.frame("Estimate" = round(unname(summary(mr.coda[, "theta"])$statistics[1]), 3),
  #                          "SD" = round(unname(summary(mr.coda[, "theta"])$statistics[2]), 3),
  #                          "2.5% quantile" = round(unname(summary(mr.coda[, "theta"])$quantiles[1]), 3),
  #                          "97.5% quantile" = round(unname(summary(mr.coda[, "theta"])$quantiles[5]), 3),
  #                          "Rhat" = round(unname(gelman.diag(mr.coda)$psrf[1]), 3))
  # names(mr_estimate) = c("Estimate", "SD", "2.5% quantile", "97.5% quantile", "Rhat")
  # return(list("MR_Estimate" = mr_estimate, "MR_Coda" = mr.coda))
}
