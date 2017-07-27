data {
  int<lower=1> N;
  real X[N];
}

generated quantities {
  real beta;
  real alpha;
  real<lower=0> sigma;
  real y[N];

  beta = normal_rng(0, 10);
  alpha = normal_rng(0, 10);
  sigma = fabs(normal_rng(0, 5));

  for (n in 1:N)
    y[n] = normal_rng(X[n] * beta + alpha, sigma);
}

