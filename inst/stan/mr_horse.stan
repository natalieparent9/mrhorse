
// Run rstantools::rstan_config() and devtools::load_all() after modifying this file to trigger recompilation, sometimes install() is also required
// In src folder the stanExports_mr_horse.o file should be updated
// Note: normal() expects sd and multi_normal expects covar matrix with variance

data {
  int<lower=0> N;                 // Number of variants (rows)
  matrix[N,2] obs;                // Observed bx and by
  vector[N] sy;                   // Standard errors for by
  vector[N] sx;                   // Standard errors for bx
  real<lower=-1> fixed_tau;       // Fixed tau value
  real<lower=-1, upper=1> omega;  // correlation parameter for bx and by
}

parameters {
  real theta;                     // Causal effect of X on Y
  vector[N] alpha;                // Pleiotropic effects on the outcome
  vector[N] bx0;                  // Latent effect of the variant on the exposure
  vector<lower=0>[N] a;           // Parameter for phi[i]
  vector<lower=0>[N] b;           // Parameter for phi[i]
  vector<lower=0, upper=1>[N] r;  // Correlation between alpha and bx0
  real<lower=0> c;                // Parameter for tau
  real<lower=0> d;                // Parameter for tau
  real<lower=0> vx0;              // Variance for bx0 distribution
  real mx0;                       // Mean for bx0 distribution
}

transformed parameters {
  vector<lower=0>[N] phi = a ./ sqrt(b);      // Scaling factor for each variant based on parameters a and b
  vector[N] rho = 2 * r - 1;                  // Converts the truncated beta distribution parameter r[i] to a correlation parameter rho[i]
  real<lower=0> tau = c / sqrt(d);            // Controls the global level of shrinkage for alphas
  if (fixed_tau != -1) tau = fixed_tau;       // If default value of -1 is given, estimate tau, otherwise fix

  matrix[N,2] mu;                             // Expected means
  mu[:, 1] = bx0;                             // Expected mean of bx
  mu[:, 2] = theta * bx0 + alpha;             // Expected mean of by

  cov_matrix[2] Sigma[N];                     // An array of N 2x2 covariance matrices
  for (i in 1:N) {
    Sigma[i][1,1] = square(sx[i]);            // Variance bx
    Sigma[i][2,2] = square(sy[i]);            // Variance by
    Sigma[i][1,2] = omega * sx[i] * sy[i];    // Covariance
    Sigma[i][2,1] = Sigma[i][1,2];
  }
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

  for (i in 1:N){
    obs[i] ~ multi_normal(mu[i], Sigma[i]);   // Jointly model bx and by likelihood
  }
}
