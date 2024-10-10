
data {
  int<lower=0> N;        // Number of variants (rows)
  vector[N] by;          // Association between the variants and exposure
  vector[N] bx;          // Association between the exposure and outcome
  vector[N] sy;          // Standard errors for by
  vector[N] sx;          // Standard errors for bx
}

parameters {
  real theta;            // Causal effect of X on Y
  vector[N] alpha;       // Pleiotropic effects on the outcome
  vector[N] bx0;         // Latent effect of the variant on the exposure
  vector<lower=0>[N] a;  // Parameter for phi[i]
  vector<lower=0>[N] b;  // Parameter for phi[i]
  vector<lower=0, upper=1>[N] r;  // Correlation between alpha and bx0
  real<lower=0> c;       // Parameter for tau
  real<lower=0> d;       // Parameter for tau
  real<lower=0> vx0;     // Variance for bx0 distribution
  real mx0;              // Mean for bx0 distribution
}

transformed parameters {
  vector<lower=0>[N] phi = a ./ sqrt(b);  // Scaling factor for each variant based on parameters a and b
  vector[N] rho = 2 * r - 1;              // Converts the truncated beta distribution parameter r[i] to a correlation parameter rho[i]
  vector[N] mu = theta * bx0 + alpha;     // Computes the mean effect on Y for each variant
  real<lower=0> tau = c / sqrt(d);        // Controls the global level of shrinkage for alphas

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

  // Likelihood
  by ~ normal(mu, sy);                         // Likelihood for by (observed beta_y)
  bx ~ normal(bx0, sx);                        // Likelihood for bx (observed beta_x)
}
