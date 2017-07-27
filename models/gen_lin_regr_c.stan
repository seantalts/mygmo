data {
  int<lower=1> N;
  real X[N];
}

generated quantities {
  real beta;
  real alpha;
  real y[N];

  beta = normal_rng(0, 10);
  alpha = normal_rng(0, 10);

  for (n in 1:N)
    y[n] = normal_rng(X[n] * beta + alpha, 1.2);
}

