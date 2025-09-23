#Start Stock/Recruit analysis

#paths----
base_path<-getwd() #May have to manually change this if you don't open as folder 
dat_path<-paste0(base_path,"/Data/")
code_pat<-paste0(base_path,"/Scripts/")
out_path<-paste0(base_path,"/Output/")

#Packages----
require("RTMB")

#Load data----
#Catch/Survey data
catch_mat<-readRDS(paste0(dat_path,"catch_mat.rds"))
GN_mat<-readRDS(paste0(dat_path,"GN_mat.rds"))
offGN_mat<-readRDS(paste0(dat_path,"offGN_mat.rds"))
LW_array<-readRDS(paste0(dat_path,"LW_byage.rds"))
M<-readRDS(paste0(dat_path,"LorenzenM.rds"))
opt_vpa<-readRDS(paste0(dat_path,"opt_vpa.rds"))
mat_dat<-readRDS(file = paste0(dat_path,"mat_byage.rds"))

#Calc. SSB----
#I'm going backcast weight-at-age and maturity at age from 1999 to 1985.  Assume weight-at-age, and maturity-at-age from 1986-1998 = Mean values from 1999-2004
hist_weights<-apply(LW_array[as.character(1999:2004),,"weight"],MARGIN = 2,FUN = mean)
hist_mat<-apply(mat_dat[as.character(1999:2004),],MARGIN = 2,FUN = mean)

hist_weights<-matrix(hist_weights,nrow = length(1986:1998),ncol=dim(LW_array)[2],byrow=TRUE,dimnames=list(seq(1986,1998),dimnames(LW_array)[[2]]))
hist_mat<-matrix(hist_mat,nrow = length(1986:1998),ncol=dim(mat_dat)[2],byrow=TRUE,dimnames=list(seq(1986,1998),dimnames(mat_dat)[[2]]))

all_weights<-rbind(hist_weights,LW_array[,,"weight"])
all_mat<-rbind(hist_mat,mat_dat)

SSB_mat<-opt_vpa$N_mat*t(all_weights[,as.character(1:6)])*0.00220462*t(all_mat[,as.character(1:6)]) #1g = 0.00220462lbs

#Plot log(R/SSB) vs SSB ----
#Note need to shift this all over by a year.

R<-opt_vpa$N_mat["age1",]/1000
SSB<-apply(X = SSB_mat,MARGIN = 2,FUN = sum)/1000

#Trim R and SSB to reflect age-1 lag
R<-R[as.character(1987:2024)]
SSB<-SSB[as.character(1986:2023)]

#Remove last 4 years from analysis due to uncertainty of those estimtes
R<-R[1:(length(R)-4)]
SSB<-SSB[1:(length(SSB)-4)] 

plot(0,0,type="n",xlim=c(0,2200),ylim=c(-2,3),las=2,xaxt="n",xlab="SSB (x1,000 lbs)",ylab="log(N_age1/SSB)")
  axis(1)
  for(i in 2:length(SSB)){
    arrows(x0 = SSB[i-1],y0 = log(R/SSB)[i-1],x1 = SSB[i],y1 = log(R/SSB)[i],length = 0.2)
    }
  points(SSB,log(R/SSB),pch=21,col="black",bg="darkgrey")
text(SSB,log(R/SSB),labels = names(SSB),adj = c(0,1))  

#Plot S/R relationship
par(mar=c(4.5,5,2,2))
plot(SSB,R,type="p",pch=21,col="black",bg="darkgrey",cex=1.5,las=2,xaxt="n",ylab="",xlab="SSB (x1,000 lbs)",xlim=c(0,2250))
  axis(1)
  mtext(text = "# Age-1 (x1,000)",side = 2,line = 4)  
text(SSB,R,labels = names(SSB),adj = c(0,1),pos = 4)

#Zoom in 
par(mar=c(4.5,5,2,2))
plot(SSB,R,type="p",pch=21,col="black",bg="darkgrey",cex=1.5,las=2,xaxt="n",ylab="",xlab="SSB (x1,000 lbs)",xlim=c(0,500),ylim=c(0,1000))
  axis(1)
  mtext(text = "# Age-1 (x1,000)",side = 2,line = 4)  
text(SSB,R,labels = names(SSB),adj = c(0,1),pos = 4)

#woof.

#Stock Recruit Modeling----
#For fun try fitting the depensatory Ricker to the full dataset
##Depensatory SR----
###Ricker denpen----
Ricker_depen<-function(params){
  getAll(dater,params)
  ln_pred_R<-2*log(S)-log(S+h)+ln_alpha-exp(ln_beta)*S
  log(R) %~% dnorm(ln_pred_R,exp(ln_sigma))
  }

params<-list(ln_alpha = 1,ln_beta = 0,ln_sigma =50, h = 1500)
dater<-list(S = SSB, R = R)

Ricker_depen(params = params)
obj<-MakeADFun(func = Ricker_depen,parameters = params)
opt<-nlminb(start = obj$par,objective =  obj$fn,gradient =  obj$gr)
opt
sdr<-sdreport(obj = obj)
sdr #standard errors on ln_alpha and h are not good

#Plot results
#Plot S/R relationship
par(mar=c(4.5,5,2,2))
  plot(SSB,R,type="p",pch=21,col="black",bg="darkgrey",cex=1.5,las=2,xaxt="n",ylab="",xlab="SSB (x1,000 lbs)",xlim=c(0,6000))
  axis(1)
  mtext(text = "# Age-1 (x1,000)",side = 2,line = 4)  
  text(SSB,R,labels = names(SSB),adj = c(0,1),pos = 4)

  est_h<-sdr$par.fixed["h"]
  est_alpha<-exp(sdr$par.fixed["ln_alpha"]+0.5*sqrt(sdr$cov.fixed["ln_sigma","ln_sigma"]))
  est_beta<-exp(sdr$par.fixed["ln_beta"])
curve((x/(x+est_h))*x*est_alpha*exp(-est_beta*x),from=c(0,2250),lty=1,lwd=2,col="red",add=TRUE)

###BH Depen----
BH_depen<-function(params){
  getAll(dater,params)
  ln_pred_R<-2*log(S)-log(S+h)+ln_alpha-log(1+exp(ln_beta)*S)
  log(R) %~% dnorm(ln_pred_R,exp(ln_sigma))
  }

params<-list(ln_alpha = 1,ln_beta = 10,ln_sigma =50, h = 1500)
dater<-list(S = SSB, R = R)

BH_depen(params = params)
obj<-MakeADFun(func = BH_depen,parameters = params)
opt<-nlminb(start = obj$par,objective =  obj$fn,gradient =  obj$gr)
opt
sdr<-sdreport(obj = obj)
sdr #standard errors on ln_alpha and h are still not good

#Plot S/R relationship still not good.
par(mar=c(4.5,5,2,2))
  plot(SSB,R,type="p",pch=21,col="black",bg="darkgrey",cex=1.5,las=2,xaxt="n",ylab="",xlab="SSB (x1,000 lbs)",xlim=c(0,2500),ylim=c(0,30000))
    axis(1)
    mtext(text = "# Age-1 (x1,000)",side = 2,line = 4)  
  text(SSB,R,labels = names(SSB),adj = c(0,1),pos = 4)

  est_h<-sdr$par.fixed["h"]
  est_alpha<-exp(sdr$par.fixed["ln_alpha"]+0.5*sqrt(sdr$cov.fixed["ln_sigma","ln_sigma"]))
  est_beta<-exp(sdr$par.fixed["ln_beta"])
curve((x/(x+est_h))*((est_alpha*x)/(1+est_beta*x)),from=c(0,2250),lty=1,lwd=2,col="red",add=TRUE)


