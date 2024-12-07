
// Run rstantools::rstan_config() and devtools::load_all() after modifying this file to trigger recompilation, sometimes install() is also required
// In src folder the stanExports_mr_horse.o file should be updated
// Note: normal() expects sd and multi_normal expects covar matrix with variance

data {
  int<lower=0> J;                 // Number of variants (rows)
  matrix[J,2] obs;                // Observed bx and by
  matrix[2,2] V[J];               // Variance covariance matrix for bx by
  real<lower=-1> fixed_tau;       // Fixed tau value
}

parameters {
  real theta;                     // Causal effect of X on Y
  vector[J] alpha;                // Pleiotropic effects on the outcome
  vector[J] bx0;                  // Latent effect of the variant on the exposure
  vector<lower=0>[J] a;           // Parameter for phi[i]
  vector<lower=0>[J] b;           // Parameter for phi[i]
  vector<lower=0, upper=1>[J] r; // Correlation between alpha and bx0
  real<lower=0> c;                // Parameter for tau
  real<lower=0> d;                // Parameter for tau
  real<lower=0> vx0;              // Variance for bx0 distribution
  real mx0;                       // Mean for bx0 distribution
}

transformed parameters {
  vector<lower=0>[J] phi = a ./ sqrt(b);      // Scaling factor for each variant based on parameters a and b
  vector[J] rho = 2 * r - 1;                  // Converts the truncated beta distribution parameter r[i] to a correlation parameter rho[i]
  real<lower=0> tau = c / sqrt(d);            // Controls the global level of shrinkage for alphas
  if (fixed_tau != -1) tau = fixed_tau;       // If default value of -1 is given, estimate tau, otherwise fix

  matrix[J,2] mu;                             // Expected means
  mu[:, 1] = bx0;                             // Expected mean of bx
  mu[:, 2] = theta * bx0 + alpha;             // Expected mean of by
}

model {
  // Priors
  theta ~ uniform(-10, 10);
  alpha ~ normal(0, tau * phi);
  bx0 ~ normal(mx0 + (sqrt(vx0) / (tau * phi)) .* rho .* alpha, sqrt((1 - rho .* rho) * vx0));
  a ~ normal(0, 1);
  b ~ gamma(0.5, 0.5);
  r ~ beta(10, 10);
  c ~ normal(0, 1);
  d ~ gamma(0.5, 0.5);
  vx0 ~ normal(0, 1);
  mx0 ~ normal(0, 1);

  for (i in 1:J){
    obs[i] ~ multi_normal(mu[i], V[i]);   // Jointly model bx and by likelihood
  }
}
