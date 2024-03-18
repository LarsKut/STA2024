data {
  int<lower=0> N; 
  vector[N] aff_i; 
  int<lower=0> observe_i[N]; 
  vector<lower=0>[N] expect_i; 
}
parameters {
  vector[N] alpha_i; 
  real beta; 
}

transformed parameters {
  vector[N] log_theta;
  log_theta = alpha_i + beta*aff_i;
}

model {
  for(i in 1:N){
    alpha_i[i] ~ normal(0,1);
  }
  beta ~ normal(0,1);
  observe_i ~ poisson_log(log_theta + log(expect_i));
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = poisson_log_lpmf(observe_i[i] | log_theta[i] + log(expect_i));
  }
}
