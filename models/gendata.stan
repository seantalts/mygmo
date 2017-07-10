parameters {
  vector<lower=0>[10] tau;
}

model {
  tau ~ cauchy(0, 2.5);
}
