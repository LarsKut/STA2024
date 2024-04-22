data {
  int<lower=1> J;         //Number of countries
  int<lower=0> T;         //Number of years
  matrix[T, J] gdp;
  matrix[T, J] urban;
  matrix[T, J] y;         //Response variable: renewable_energy
  int<lower=1> P;         //Number of years to forecast
}

parameters {
  real beta_gdp;
  real beta_urban;
  real<lower=0> sigma; 
  real<lower=-1, upper=1> rho;
}


model {
  //Priors
  beta_urban ~ normal(0,1);
  beta_gdp ~ normal(0,1);
  rho ~ uniform(-1, 1);
  sigma ~ normal(0,1);
  
  //Model response for each time point and country
  for(j in 1:J) {
    y[1,j] ~ normal(0, sigma/sqrt((1-rho^2)));
    for (t in 2:T) {
    y[t,j] ~ normal(gdp[t-1,j]*beta_gdp + urban[t-1,j]*beta_urban + rho*y[t-1,j], sigma);
    }
  }
}


generated quantities {
  matrix[T, J] y_rep;  // Replicated data for observed years
  matrix[P, J] y_pred;  // Predicted data for future years
  matrix[P, J] urban_pred; 
  matrix[P, J] gdp_pred;
  for (j in 1:J) {
  
    y_rep[1,j] = normal_rng(y[1,j], sigma / sqrt(1 - rho^2));
    for (t in 2:T) {
      y_rep[t,j] = normal_rng(gdp[t-1,j] * beta_gdp + urban[t-1,j] * beta_urban + rho * y_rep[t-1,j], sigma);
    }
    
    //Simulate the covariates using AR(1) model
    real sigma_x = 1;
    real rho_x = 0.5;
    gdp_pred[1,j] = normal_rng(gdp[T,j] * rho_x, sigma_x);
    urban_pred[1,j] = normal_rng(urban[T,j] * rho_x, sigma_x);
    
    y_pred[1,j] = normal_rng(rho * y[T,j], sigma);
    for (t in 2:P) {
      gdp_pred[t,j] = normal_rng(gdp_pred[t-1,j] * rho_x,sigma_x);
      urban_pred[t,j] = normal_rng(urban_pred[t-1,j] * rho_x, sigma_x);
      y_pred[t,j] = normal_rng(gdp_pred[t-1,j]*beta_gdp + urban_pred[t-1,j]*beta_urban + rho * y_pred[t-1,j], sigma);
    }
  }
}


