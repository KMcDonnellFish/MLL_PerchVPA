####VonB and weight/length model

VonBFUN<-function(params){
  getAll(params,vonB_dater_singleyear)

  L_pred<-Linf*(1-exp(-k*(age+0.4164-t0))) #Note this offset is to account for the ~5 months from beginning of year to time of capture during fall survey (~5 months from hatch).
  #Likelihood
  length %~% dnorm(mean = L_pred, sd = L_pred_sd)

  }