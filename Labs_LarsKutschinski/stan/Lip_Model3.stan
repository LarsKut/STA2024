data {
  int<lower=0> N; 
  vector[N] aff_i; 
  int<lower=0> observe_i[N]; 
  vector<lower=0>[N] expect_i; 
}
parameters {
  real mu_alpha; 
  real<lower=0> sigma_alpha; 
  vector[N] alpha_i; 
  real beta; 
}

transformed parameters {
  vector[N] log_theta;
  log_theta = alpha_i + beta*aff_i;
}

model {
  mu_alpha ~ normal(0,1);
  sigma_alpha ~ normal(0,1);
  alpha_i ~ normal(mu_alpha, sigma_alpha); 
  observe_i ~ poisson_log(log_theta + log(expect_i));
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = poisson_log_lpmf(observe_i[i] | log_theta[i] + log(expect_i));
  }
}
