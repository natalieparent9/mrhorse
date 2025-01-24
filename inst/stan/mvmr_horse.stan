
// After this filed is modified, it needs to be recompiled.
// Triggering recompilation can be inconsistent, first try Run rstantools::rstan_config() and devtools::load_all(), can try pkgbuild::compile_dll()
// In src folder the stanExports_mvmr_horse.o / .h file should be updated

data {
    int<lower=1> J;                    // Number of variants (rows)
    int<lower=1> K;                    // Number of exposures
    matrix[J,K+1] obs;                 // Observed bx and by
    array[J] matrix[K+1, K+1] V;       // Variance covariance matrix for bxs & by
    matrix[K, K] R;                    // Prior precision matrix for mx, identity
    real<lower=-1> fixed_tau;          // Fixed tau value
}

parameters {
    matrix[J, K] bx0;                  // Latent effect of the variants on the exposures
    vector[K] theta;                   // Causal effect of X[k] on Y
    vector[J] alpha;                   // Pleiotropic effects on the outcome
    vector<lower=0, upper=1>[J] r;     // Will be scaled to correlation between alpha and bx0 rho
    vector<lower=0>[J] a;              // Parameter for phi[i]
    vector<lower=0>[J] b;              // Parameter for phi[i]
    real<lower=0> c;                   // Parameter for tau
    real<lower=0> d;                   // Parameter for tau
    vector[K] mx;                      // Mean for bx0 distribution
    vector<lower=0>[K] vx0;            // Variance for bx0 distribution
}

transformed parameters {
    vector[J] phi = a ./ sqrt(b);                              // Scaling factor for each variant based on parameters a and b
    vector[J] rho = 2 * r - 1;                                 // Converts the truncated beta distribution parameter r[i] to a correlation parameter rho[i]
    vector[J] kappa = (rho^2 ./ (1 + K*rho^2));                // Used to adjust bx0
    matrix[K, K] A = diag_matrix(1.0 ./ vx0);                  // Diagonal precision matrix for bx0, diagonal elements are 1/vx0, off-diagonal are 0
    matrix[K, K] B = (1.0 ./ sqrt(vx0)) * (1.0 ./ sqrt(vx0))'; // Matrix B for covariance
    real<lower=0> tau = c / sqrt(d);                           // Controls the global level of shrinkage for alphas
    if (fixed_tau != -1) tau = fixed_tau;                      // If default value of -1 is given, estimate tau, otherwise fix
    matrix[K, K] precision_matrix[J];
    matrix[J,K+1] mu;                                          // Expected means
    for (i in 1:J) {
      precision_matrix[i] = A - kappa[i] * B;
      for (k in 1:K) mu[i,k] = bx0[i,k];                      // Expected mean of bxs
      mu[i, K+1] = dot_product(bx0[i], theta) + alpha[i];     // Expected mean of by
    }
}

model {
    for (i in 1:J) {
        // Likelihood
        obs[i] ~ multi_normal(mu[i], V[i]);

        // Priors
        alpha[i] ~ normal(0, tau * phi[i]);
        bx0[i] ~ multi_normal(mx + sqrt(vx0) .* (rho[i] * alpha[i] / (phi[i] * tau)), inverse(precision_matrix[i]));
        r[i] ~ beta(10, 10);
        a[i] ~ normal(0, 1);
        b[i] ~ gamma(0.5, 0.5);
    }
    // Priors
    c ~ normal(0, 1);
    d ~ gamma(0.5, 0.5);
    mx ~ multi_normal(rep_vector(0, K), R);

    // Parameters estimated for each exposure
    for (k in 1:K) {
        vx0[k] ~ normal(0, 1);
        theta[k] ~ uniform(-10, 10);
    }
}

