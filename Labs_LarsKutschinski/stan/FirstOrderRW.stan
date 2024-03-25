data {
  int<lower=0> N; // number of observations
  int<lower=0> T; //number of years
  int<lower=0> mid_year; // mid-year of study
  vector[N] y; //log ratio
  int<lower=0> P;
  vector[N] se; // standard error around observations
  vector[T] years; // unique years of study
  int<lower=0> year_i[N]; // year index of observations
  
}

parameters {
  real<lower=0> sigma_mu; 
  vector[T] mu; 
  
}

model {

  y ~ normal(mu[year_i], se);
  sigma_mu ~ normal(0,1);
  mu[1]~normal(0,sigma_mu);
  mu[2:T]~normal(mu[1:(T-1)],sigma_mu);
}

generated quantities {
  vector[P] projected_mu;
  projected_mu[1] = normal_rng(mu[T], sigma_mu);
  for (i in 2:P) {
    projected_mu[i] = normal_rng(projected_mu[i-1], sigma_mu);
  }
}