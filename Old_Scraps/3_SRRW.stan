data{
  int n;
  real S[n];
  real log_R[n];
  }

parameters{
  real log_alpha;
  real<lower=0> beta;
  vector[n] phi;
  real<lower=0> phi_sd;
  real<lower=0> eps_sd;
  #real w0;
  }

model{
  vector[n] w;
  vector[n] log_pred_R;
  
  #Priors
  log_alpha ~ normal(0,10);
  #beta ~ uniform(0,1000); #let stan assign wide uniform
  eps_sd ~ gamma(2,0.1);#normal(0,10);
  #w0 ~ normal(0,20);
  phi_sd ~ gamma(2,0.1);#normal(0,10);
  
  phi~normal(0,phi_sd);
  
  #Model
  #Initiate Predictions
  #w[1] = w0;
  #log_pred_R[1] = log(S[1]) + log_alpha - beta*S[1] + w[1];
  
  #Generate RW predictions
  for(i in 1:n){
    if(i==1){
      w[1] = phi[1];
      } else {
        w[i] = w[(i-1)] + phi[i];
        }
    log_pred_R[i] = log(S[i]) + log_alpha - beta*S[i] + w[i];
    }
  
  #Likelihood Eval
  log_R ~ normal(log_pred_R,eps_sd);
  }
generated quantities{
  vector[n] w;
  vector[n] log_pred_R;
  vector[n] log_alpha_w;

  #Generate RW predictions
  for(i in 1:n){
    if(i==1){
      w[1] = phi[1];
      } else {
        w[i] = w[(i-1)] + phi[i];
        }
    log_pred_R[i] = log(S[i]) + log_alpha - beta*S[i] + w[i];
    }
  
  log_alpha_w = log_alpha + w;
  
  }
