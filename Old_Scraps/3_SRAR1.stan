data{
  int n;
  real S[n];
  real log_R[n];
  }

parameters{
  real log_alpha;
  real<lower=0> beta;
  real<lower=(-1),upper=1> phi;
  real<lower=0> eps_sd;
  #real<lower=0> log_mean_R0;
  #real<lower=0> sd_R0;
  real log_pred_R0;
  }

model{
  vector[n] log_pred_R;
  vector[(n-1)] w;
  
  #Priors
  log_alpha ~ uniform(-20,20);#normal(0,10);
  #beta ~ uniform(0,1000); #let stan assign wide uniform
  phi ~ uniform(-1,1);
  eps_sd ~ normal(0,10);
  #log_mean_R0 ~ normal(0,50); #let stan assign wide uniform
  #sd_R0 ~ normal(0,10);
  #log_pred_R0~normal(log_mean_R0,sd_R0);
  log_pred_R0~normal(12.64562,12.64562*0.4); #Informative prior based on mean value from first 5 years of dataset.
  
  #Model
  #Initiate Predictions
  log_pred_R[1] = log_pred_R0;
  
  #Generate AR1 predictions
  for(i in 2:n){
    w[(i-1)] = log_pred_R[(i-1)] - (log(S[(i-1)]) + log_alpha - beta*S[(i-1)]);
    log_pred_R[i] = log(S[i]) + log_alpha - beta*S[i] + phi*w[(i-1)];
    }
  
  #Likelihood Eval
  log_R ~ normal(log_pred_R,eps_sd);
  }

