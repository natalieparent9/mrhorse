data {
    int<lower=1> N;                    // Number of observations
    int<lower=1> K;                    // Number of variants
    matrix[N, K] bx;                   // Matrix of observed effect estimates on x
    vector[N] by;                      // Vector of observed effect estimates on y
    vector[N] sy;                      // Standard errors for by
    matrix[K, K*N] Tx;                 // Covariance matrix for bx
    matrix[K, K] R;                    // Prior precision matrix for mx, diagonal of 1s
}

parameters {
    matrix[N, K] bx0;                  // Latent matrix of effect estimates on x
    vector[K] theta;                   // Effect estimates on y
    vector[N] alpha;                   // Random intercepts
    vector<lower=0, upper=1>[N] r;     // Beta distribution parameter for rho
    vector<lower=0>[N] a;              // Parameter for phi
    vector<lower=0>[N] b;              // Parameter for phi
    vector[K] mx;                      // Prior mean for bx0
    vector<lower=0>[K] vx0;            // Variances for bx0
    real<lower=0> c;                   // Parameter for tau
    real<lower=0> d;                   // Parameter for tau
}

transformed parameters {
    real tau = c / sqrt(d);                  // Tau is defined as c divided by sqrt(d)
    vector[N] phi = a ./ sqrt(b);            // Scaling factor for each SNP based on parameters a and b, ./ is element wise division
    vector[N] rho = 2 * r - 1;               // Converts the truncated beta distribution parameter r[i] to a correlation parameter rho[i]
    vector[N] kappa = (rho^2 ./ (1 + K*rho^2));   // Used to adjust bx0
    matrix[K, K] A = diag_matrix(1.0 ./ vx0);  // Diagonal precision matrix for bx0, diagonal elements are 1/vx0, off-diagonal are 0
    matrix[K, K] B = (1.0 ./ sqrt(vx0)) * (1.0 ./ sqrt(vx0))'; // Matrix B for covariance
    matrix[K, K] covariance_matrix[N];
    // Convert precision matrix to covariance
    for (i in 1:N) {
      matrix[K, K] precision_matrix = A - kappa[i] * B; // covariance matrix for bx0
      covariance_matrix[i] = inverse(precision_matrix);
    }
}


model {
    for (i in 1:N) {
        // Likelihood
        by[i] ~ normal(dot_product(bx0[i], theta) + alpha[i], sy[i]);
        bx[i] ~ multi_normal(bx0[i], Tx[, (i-1)*K + 1:i*K]);

        // Priors
        alpha[i] ~ normal(0, tau * phi[i]);
        bx0[i] ~ multi_normal(mx + sqrt(vx0) .* (rho[i] * alpha[i] / (phi[i] * tau)), covariance_matrix[i]);
        r[i] ~ beta(10, 10);
        a[i] ~ normal(0, 1);
        b[i] ~ gamma(0.5, 0.5);
    }

    // Priors
    c ~ normal(0, 1);
    d ~ gamma(0.5, 0.5);

    mx ~ multi_normal(rep_vector(0, K), R);

    for (k in 1:K) {
        vx0[k] ~ normal(0, 1);
        theta[k] ~ uniform(-10, 10);
    }
}

