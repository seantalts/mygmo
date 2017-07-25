data {
  int<lower=1> N;
  vector[N] X;
  vector[N] y;
}

parameters {
  real beta;
  real alpha;
  real<lower=0> sigma;
}

model {
  beta ~ normal(0, 10);
  alpha ~ normal(0, 10);
  sigma ~ normal(0, 5);

  y ~ normal(X * beta + alpha, sigma);
}
