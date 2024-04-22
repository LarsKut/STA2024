data {
  int<lower=1> J;                    // Number of countries
  int<lower=0> T;                    //number of years
  matrix[T, J] gdp;
  matrix[T, J] urban;
  matrix[T, J] y;                       // Response variable: renewable_energy
  int<lower=1> P;
}

parameters {
  real beta_gdp;
  real beta_urban;
  real<lower=0> sigma; 
  real<lower=-1,upper=1> rho_1;
  real<lower=-1,upper=1> rho_2;
}



model {
  beta_urban ~ normal(0,1);
  beta_gdp ~ normal(0,1);
  rho_1 ~ uniform(-0.5,0.5);
  rho_2 ~ uniform(-0.5,0.5);
  sigma ~ normal(0,1);
  
  for(j in 1:J) {
    y[1,j] ~ normal(y[1,j], sigma/sqrt((1-rho_1^2)));
    y[2,j] ~ normal(y[2,j], sigma/sqrt((1-rho_2^2)));
    for (t in 3:T) {
    y[t,j] ~ normal(gdp[t-1,j]*beta_gdp + urban[t-1,j]*beta_urban + rho_1*y[t-1,j] + rho_2*y[t-2,j], sigma);
    }
  }
}


generated quantities {
  matrix[T, J] y_rep;
  matrix[P, J] y_pred;  // Replicated data for posterior checks
 
  for (j in 1:J) {
    
    y_rep[1,j] = normal_rng(y[1,j], sigma / sqrt(1 - rho_1^2));
    y_rep[2,j] = normal_rng(y[2,j], sigma / sqrt(1 - rho_2^2));
    // Simulate subsequent observed years
    for (t in 3:T) {
      y_rep[t,j] = normal_rng(gdp[t-1,j] * beta_gdp + urban[t-1,j] * beta_urban + rho_1 * y_rep[t-1,j] + rho_2 * y_rep[t-2,j], sigma);
    }
    
    // Simulate the first observation
    y_pred[1,j] = normal_rng(rho_1 * y[T,j] + rho_2 * y[T-1,j], sigma);
    y_pred[2,j] = normal_rng(rho_1 * y_pred[1,j] + rho_2 * y[T,j], sigma);

    // Simulate subsequent observations
    for (t in 3:P) {
     
      y_pred[t,j] = normal_rng(rho_1*y_pred[t-1,j] + rho_2*y_pred[t-2,j], sigma);
    }
  }
}
