data{
  int n;
  vector[n] S;
  vector[n] log_R;
  }

parameters{
//  vector[(n-1)] phi;
//  real<lower=0> phi_sd;
  real log_alpha;
  real<lower=0> beta;
  real<lower=0,upper=500000> h;
  real<lower=0> eps_sd;

  }

model{
  vector[n] log_pred_R;

  #Priors
  eps_sd ~ gamma(2,1/0.1);#normal(0,10);
  
  #Calculate expected vals
  log_pred_R = 2*log(S)-log(S+h)+log_alpha-log(1+beta*S);
  
  #Likelihood Eval
  log_R ~ normal(log_pred_R,eps_sd);
  }