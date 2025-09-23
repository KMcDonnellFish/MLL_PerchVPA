data{
  int n;
  real S[n];
  real log_R[n];
  real sbr0[n];
  }

parameters{
//  vector[(n-1)] phi;
//  real<lower=0> phi_sd;
  real log_recK;
  real log_R0;
  real<lower=0> eps_sd;

  }

model{
  //vector[n] w;
  vector[n] log_pred_R;
//  vector[n] log_CR;
  
  #Priors
//  phi_sd ~ gamma(2,1/0.1);#normal(0,10);
//  phi ~ normal(0,phi_sd);
  log_recK ~ normal(2.287,0.893);
  //log_R0 ~ gamma(2,1/0.1);
  eps_sd ~ gamma(2,1/0.1);#normal(0,10);
  
  #Model
  #Generate RW predictions
  for(i in 1:n){
//    if(i==1){
//      log_CR[1] = log_recK;
//      } else {
//        log_CR[i] = log_CR[(i-1)] + phi[(i-1)];
//        }
    
//    log_CR = log_recK;
    
    #log scale Ricker    
    log_pred_R[i] = log(S[i])-log(sbr0[i])+log_recK*(1-(S[i]/(exp(log_R0)*sbr0[i])));
    }
  
  #Likelihood Eval
  log_R ~ normal(log_pred_R,eps_sd);
  }