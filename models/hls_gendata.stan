data {
  int<lower=0> N;              // num individuals
  int<lower=1> K;              // num ind predictors
  int<lower=1> J;              // num groups
  int<lower=1> L;              // num group predictors
  int<lower=1,upper=J> jj[N];  // group for individual
  matrix[N, K] x;               // individual predictors
  matrix[J, L] u;              // group predictors
}

generated quantities {
  matrix[K, J] z;
  cholesky_factor_corr[K] L_Omega;
  matrix[K, K] Omega;             // prior correlation
  vector<lower=0>[K] tau;      // prior scale
  matrix[L, K] gamma;           // group coeffs
  matrix[J, K] beta;           // indiv coeffs by group
  real<lower=0> sigma;         // prediction error scale
  vector[N] y;                 // outcomes

  L_Omega = lkj_corr_cholesky_rng(K, 4);

  for (k in 1:K) {
    tau[k] = -1;
    while (tau[k] <= 0)
      tau[k] = normal_rng(0, 5);
    for (j in 1:J)
      z[k, j] = normal_rng(0, 1);
  }

  for (i in 1:L)
    for (j in 1:K)
      gamma[i, j] = normal_rng(0, 5);

  sigma = -1;
  while (sigma <= 0)
    sigma = normal_rng(0, 5);

  beta = u * gamma + (diag_pre_multiply(tau, L_Omega) * z)';
  
  for (n in 1:N)
    y[n] = normal_rng(dot_product(beta[jj[n]], x[n]), sigma);
  
  Omega = L_Omega' * L_Omega;
}
