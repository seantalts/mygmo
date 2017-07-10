data {
  int<lower=0> N;              // num individuals
  int<lower=1> K;              // num ind predictors
  int<lower=1> J;              // num groups
  int<lower=1> L;              // num group predictors
  int<lower=1,upper=J> jj[N];  // group for individual
}

generated quantities {
  corr_matrix[K] Omega;        // prior correlation
  vector<lower=0>[K] tau;      // prior scale
  matrix[L, K] gamma;           // group coeffs
  vector[K] beta[J];           // indiv coeffs by group
  real<lower=0> sigma;         // prediction error scale
  matrix[N, K] x;               // individual predictors
  row_vector[L] u[J];          // group predictors
  vector[N] y;                 // outcomes

  Omega = lkj_corr_rng(K, 2);
  for (n in 1:K) {
    tau[n] = fabs(cauchy_rng(0, 2.5));
  }
  for (i in 1:L) {
    for (j in 1:K) {  
      gamma[i, j] = normal_rng(0, 5);
    }
  }
  sigma = fabs(cauchy_rng(0, 2.5));
  for (j in 1:J)
    for (k in 1:K)
      beta[j, k] = normal_rng(0, 5);


  for (n in 1:N)
    for (k in 1:K)
      x[n, k] = normal_rng(0, 10);
      
  for (j in 1:J)
    for (l in 1:L)
      u[j, l] = normal_rng(0, 5);
  {
    row_vector[K] u_gamma[J];
    for (j in 1:J) {
      u_gamma[j] = u[j] * gamma;
      beta[j] = multi_normal_rng(to_vector(u_gamma[j]), quad_form_diag(Omega, tau));
    }
  }
  for (n in 1:N)
    y[n] = normal_rng(x[n] * beta[jj[n]], sigma);
}
