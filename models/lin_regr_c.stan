data {
  int<lower=1> N;
  vector[N] X;
  vector[N] y;
}

parameters {
  real beta;
  real alpha;
}

model {
  beta ~ normal(0, 10);
  alpha ~ normal(0, 10);

  y ~ normal(X * beta + alpha, 1.2);
}
