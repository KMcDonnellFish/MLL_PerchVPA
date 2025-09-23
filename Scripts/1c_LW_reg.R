LW_reg<-function(params){
  getAll(vonB_dater_singleyear,params)
  
  #log_a %~% dnorm(mean = 0,sd = 2)
  #b %~% dnorm(mean = 0,sd = 10)
  #w_sd %~% dgamma(shape = 0.001,scale = 0.001)
  
  pred_weight<-exp(log_a+b*log(length))
  
  weight %~% dnorm(mean = pred_weight,sd = w_sd)
  
  }