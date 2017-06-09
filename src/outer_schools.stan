data {
  int<lower=0> J; // number of schools
  real y[J]; // estimated treatment effects
  real<lower=0> sigma[J]; // s.e. of effect estimates
  real eta[J]; // school-level errors - guess of local parameters
}
parameters {
  real mu; // population mean
  real<lower=0> tau; // population sd
}
transformed parameters {
  real theta[J]; // school effects
  for (j in 1:J)
    theta[j] = mu + tau * eta[j];
}
model {
  eta ~ normal(0, 1);
  y ~ normal(theta, sigma);
}
