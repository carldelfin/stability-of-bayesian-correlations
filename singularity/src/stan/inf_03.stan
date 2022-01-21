data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}
parameters {
  real<lower=-1, upper=1> beta;
  real<lower=0> sigma;
}
transformed parameters {
  real<lower=0> beta_transformed;
  beta_transformed = (beta + 1) / 2;
}
model {
  beta  ~ normal(0, 2.5);
  sigma ~ cauchy(0, 1);
  target += beta_lpdf(beta_transformed | 8.8, 5.2);
  y ~ normal(0 + beta * x, sigma);
}

