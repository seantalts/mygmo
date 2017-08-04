data {
  int<lower=0> J;          // number of schools
  int<lower=1> K;          // number of hyperparameters
  real<lower=0> sigma[J];  // std err of effect estimate (school j)
}

generated quantities {
  real mu;
  real<lower=0> tau;
  real theta[J];
  
  // data to generate:
  real y[J];               // estimated treatment effect (school j)
  
  mu = normal_rng(0, 5);
  tau = fabs(normal_rng(0, 4));
  //tau = -1;
  //while (tau <= 0)
  //  tau = normal_rng(0, 5);

  for (j in 1:J) {
    theta[j] = normal_rng(mu, tau);
    y[j] = normal_rng(theta[j], sigma[j]);
  }
}
