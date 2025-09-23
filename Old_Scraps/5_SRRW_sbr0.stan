data{
  int n;
  real S[n];
  real log_R[n];
  real sbr0[n];
  }

parameters{
  vector[n-1] phi;
  real<lower=0> phi_sd;
  real log_recK;
  real log_R0;
  real<lower=0> eps_sd;

  }

model{
  vector[n] w;
  vector[n] CR;
  vector[n] alpha;
  vector[n] beta;
  vector[n] log_pred_R;
  
  #Priors
  phi_sd ~ gamma(2,1/0.1);#normal(0,10);
  phi ~ normal(0,phi_sd);
  log_recK ~ normal(0,10);
  log_R0 ~ uniform(0,100);
  eps_sd ~ gamma(2,1/0.1);#normal(0,10);
  
  #Model
  #Generate RW predictions
  for(i in 1:n){
    if(i==1){
      w[1] = log_recK;
      } else {
        w[i] = w[(i-1)] + phi[(i-1)];
        }
    
    #Create easier intermediates
    CR[i] = exp(w[i]);
    alpha[i] = CR[i]/sbr0[i];
    beta[i] = log(CR[i])/(exp(log_R0)*sbr0[i]);
    #log scale Ricker    
    log_pred_R[i] = log(alpha[i]*S[i]*exp(-beta[i]*S[i]));
    }
  
  #Likelihood Eval
  log_R ~ normal(log_pred_R,eps_sd);
  }
// generated quantities{
//   vector[n] w;
//   vector[n] CR;
//   vector[n] alpha;
//   vector[n] beta;
//   vector[n] log_pred_R;
//   
//   #Model
//   #Generate RW predictions
//   for(i in 1:n){
//     if(i==1){
//       w[1] = phi[1];
//       } else {
//         w[i] = w[(i-1)] + phi[i];
//         }
//     
//     #Create easier intermediates
//     CR[i] = exp(log_recK+w[i]);    
//     alpha[i] = CR[i]/sbr0[i];
//     beta[i] = log(CR[i])/(exp(log_R0)*sbr0[i]);
//     #log scale Ricker    
//     log_pred_R[i] = log(alpha[i]*S[i]*exp(-beta[i]*S[i]));
//     }
//   }