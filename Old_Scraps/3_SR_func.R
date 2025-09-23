#Build a RTMB function to fit S/R models

KMc_SR<-function(params){
  getAll(dater,params,model_opt)
  
  S<-SR_dater[,"SSB"]
  
  if(mod_type==1){ #BH - params = alpha, beta, sigma
    est_log_recruit<-log((exp(log_alpha)*S)/(exp(log_beta)+S))
    }
  if(mod_type==2){ #Ricker - params = alpha, log_beta, sigma
    est_log_recruit<-log(alpha*S*exp(-exp(log_beta)*S))
    }
  if(mod_type==3){
    #est_log_recruit<-log(exp(log_alpha)*(S^gamma)/(1+exp(log_beta)*(S^gamma)))
    est_log_recruit<-log((S/(S+exp(log_gamma)))*((exp(log_alpha)*S)/(1+exp(log_beta)*S)))
    }
  if(mod_type==4){
    est_log_recruit<-log((S/(S+exp(log_gamma)))*(exp(log_alpha)*S*exp(-exp(log_beta)*S)))
    }

  #Evaluate likelihood
  log(SR_dater[,"recruit"])%~%dnorm(mean = est_log_recruit,sd = sigma)
  
  }

KMc_SR_AR<-function(params){
  getAll(dater,params)
  
  S<-SR_dater[,"SSB"]
  R<-SR_dater[,"recruit"]
  
  R_pred<-R[1]
  
  for(i in 2:length(R)){
    temp_R<-exp(log_alpha)*S[i]*exp(-exp(log_beta)*S[i]) + exp(log_gamma)*R[i-1]
    R_pred<-c(R_pred,temp_R)
    }
  
  log(R[2:length(R)])%~%dnorm(mean = log(R_pred[2:length(R)]),sd = sigma)
  
  }

KMc_SR_AR_bayes<-function(params){
  getAll(dater,params)
  
  S<-SR_dater[,"SSB"]
  R<-SR_dater[,"recruit"]
  
  R_pred<-R[1]
  
  for(i in 2:length(R)){
    temp_R<-exp(log_alpha)*S[i]*exp(-exp(log_beta)*S[i]) + exp(log_gamma)*R[i-1]
    R_pred<-c(R_pred,temp_R)
    }
  
  log(R[2:length(R)])%~%dnorm(mean = log(R_pred[2:length(R)]),sd = sigma)
  
  }

  
