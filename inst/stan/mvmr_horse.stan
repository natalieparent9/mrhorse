
// Run rstantools::rstan_config() and devtools::load_all() after modifying this file to trigger recompilation
// In src folder the stanExports_mvmr_horse.o file should be updated

data {
    int<lower=1> N;                    // Number of variants (rows)
    int<lower=1> K;                    // Number of exposures
    matrix[N,K+1] obs;                 // Observed bx and by
    matrix[N,K] sx;                      // Standard errors for bx
    vector[N] sy;                      // Standard errors for by
    // matrix[K, K*N] Tx;                 // Covariance matrix for bx
    matrix[K, K] R;                    // Prior precision matrix for mx, diagonal of 1s
    real<lower=-1> fixed_tau;          // Fixed tau value
    real<lower=-1, upper=1> omega;     // Correlation parameter for bx and by
}

parameters {
    matrix[N, K] bx0;                  // Latent effect of the variants on the exposures
    vector[K] theta;                   // Causal effect of X[k] on Y
    vector[N] alpha;                   // Pleiotropic effects on the outcome
    vector<lower=0, upper=1>[N] r;     // Correlation between alpha and bx0
    vector<lower=0>[N] a;              // Parameter for phi[i]
    vector<lower=0>[N] b;              // Parameter for phi[i]
    real<lower=0> c;                   // Parameter for tau
    real<lower=0> d;                   // Parameter for tau
    vector[K] mx;                      // Mean for bx0 distribution
    vector<lower=0>[K] vx0;            // Variance for bx0 distribution
}

transformed parameters {
    vector[N] phi = a ./ sqrt(b);                              // Scaling factor for each variant based on parameters a and b
    vector[N] rho = 2 * r - 1;                                 // Converts the truncated beta distribution parameter r[i] to a correlation parameter rho[i]
    vector[N] kappa = (rho^2 ./ (1 + K*rho^2));                // Used to adjust bx0
    matrix[K, K] A = diag_matrix(1.0 ./ vx0);                  // Diagonal precision matrix for bx0, diagonal elements are 1/vx0, off-diagonal are 0
    matrix[K, K] B = (1.0 ./ sqrt(vx0)) * (1.0 ./ sqrt(vx0))'; // Matrix B for covariance
    real<lower=0> tau = c / sqrt(d);                           // Controls the global level of shrinkage for alphas
    if (fixed_tau != -1) tau = fixed_tau;                      // If default value of -1 is given, estimate tau, otherwise fix
    matrix[K, K] precision_matrix[N];
    for (i in 1:N) {
      precision_matrix[i] = A - kappa[i] * B;
    }

    matrix[N,K+1] mu;                                      // Expected means
    cov_matrix[K+1] Sigma[N];                              // An array of N K+1xK+1 covariance matrices
    for (i in 1:N) {
      Sigma[i] = rep_matrix(0, K+1, K+1);                  // First fill with zeroes
      for (k in 1:K) {
        mu[i,k] = bx0[i,k];                                // Expected mean of bxs
        Sigma[i][k,k] = square(sx[i,k]);                   // Variance of bx[k]
        Sigma[i][k,K+1] = omega * sx[i,k] * sy[i];         // Covariance
        Sigma[i][K+1,k] = Sigma[i][k,K+1];
      }
      Sigma[i][K+1,K+1] = square(sy[i]);                   // Last diagonal element, var of by
      mu[i, K+1] = dot_product(bx0[i], theta) + alpha[i];  // Expected mean of by
    }
}

model {
    for (i in 1:N) {
        // Likelihood
        obs[i] ~ multi_normal(mu[i], Sigma[i]);

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

