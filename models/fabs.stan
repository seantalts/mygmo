generated quantities {
  real fabs_normal = fabs(normal_rng(0, 1));
  real while_normal = normal_rng(0, 1);

  while (while_normal <= 0)
    while_normal = normal_rng(0, 1);
}
