data {
  int<lower=0> N;              // num individuals
  int<lower=1> K;              // num ind predictors
  int<lower=1> J;              // num groups
  int<lower=1> L;              // num group predictors
  int<lower=1,upper=J> jj[N];  // group for individual
  matrix[N, K] x;              // individual predictors
  matrix[J, L] u;              // group predictors
  vector[N] y;                 // outcomes
}
parameters {
  matrix[K, J] z;
  cholesky_factor_corr[K] L_Omega;
  vector<lower=0>[K] tau;      // prior scale
  matrix[L, K] gamma;           // group coeffs
  real<lower=0> sigma;         // prediction error scale
}
transformed parameters {
  matrix[J, K] beta;
  beta = u * gamma + (diag_pre_multiply(tau, L_Omega) * z)';
}
model {
  L_Omega ~ lkj_corr_cholesky(4);
  tau ~ normal(0, 5);
  to_vector(z) ~ normal(0, 1);
  to_vector(gamma) ~ normal(0, 5);
  sigma ~ normal(0, 5);
  y ~ normal(rows_dot_product(beta[jj] , x), sigma);
  // for logistic:
  // bernoulli_logit_regression(beta[jj], x)
}

generated quantities {
  matrix[K, K] Omega;
  Omega = L_Omega' * L_Omega;
}
