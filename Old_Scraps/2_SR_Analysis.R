
#Set Paths
dat_path<-"C:/Users/EU01237640/OneDrive - State of Minnesota - MN365/Mille_Lacs/YEP/PerchVPA/"

#I only have length/weight data from 1999 to present.  THey used scales prior to 1999 and I do not have those data.

#Load data
vpa_dat<-readRDS(file = paste0(dat_path,"opt_vpa.rds"))
LW_dat<-readRDS(file = paste0(dat_path,"LW_byage.rds"))
mat_dat<-readRDS(file = paste0(dat_path,"mat_byage.rds"))

#Going to use the age 6 weights as the plus group - need to reevaluate if I should recalc. VonBs using 6 as plus group.

#Calculate SSB
SSB<-t(vpa_dat$N_mat[,as.numeric(colnames(vpa_dat$N_mat))>=1999])*LW_dat[,as.character(seq(1,6,1)),"weight"]*mat_dat[,as.character(seq(1,6,1))]
#Convert from grams to lbs
SSB<-SSB*0.00220462

#Quick plot to look at SSB through time
plot(as.numeric(rownames(SSB)),rowSums(SSB)/1000,type="b",pch=21,col="black",bg="darkgrey",las=2,xaxt="n",ylab="SSB (x1,000lbs)",xlab="Year")
axis(1)
#Define Recruits
recruit<-vpa_dat$N_mat[1,]
recruit<-recruit[as.numeric(names(recruit))>=1999]

#Offset SSB and recruits so they match
SR_dater<-data.frame(SSB=rowSums(SSB)[-nrow(SSB)],recruit = recruit[-1])

#SR_dater<-data.frame(SSB=colSums(vpa_dat$N_mat[3:6,])[-ncol(vpa_dat$N_mat)],recruit = vpa_dat$N_mat[1,-1])

#Plot SSB vs log(age1/SSB)
plot(SR_dater[,"SSB"],log(SR_dater[,"recruit"]/SR_dater[,"SSB"]),pch=21,col="black",bg="darkgrey",las=2,xaxt="n",ylab="log(age1/SSB)",xlab="SSB (lbs)",xlim=c(0,600000))
  axis(1)
  for(i in 1:(length(SR_dater[,"recruit"])-1)){
    arrows(x0 = SR_dater[i,"SSB"],y0 = log(SR_dater[i,"recruit"]/SR_dater[i,"SSB"]),x1 = SR_dater[i+1,"SSB"],y1 = log(SR_dater[i+1,"recruit"]/SR_dater[i+1,"SSB"]),length = 0.15)
    }
  points(SR_dater[,"SSB"],log(SR_dater[,"recruit"]/SR_dater[,"SSB"]),pch=21,col="black",bg="darkgrey")
text(SR_dater[,"SSB"],log(SR_dater[,"recruit"]/SR_dater[,"SSB"]),labels = rownames(SSB),adj = c(0,1),pos = 4,cex=0.8)

#Plot out the S/R Relationship
par(mar=c(4.5,6,2,2))
plot(SR_dater[,"SSB"]/1000,SR_dater[,"recruit"]/1000,pch=21,col="black",bg="darkgrey",las=2,cex=1.2,
  xaxt="n",ylab="",xlab="SSB (x1,000 lbs)",xlim=c(0,600),ylim=c(0,1000))
  axis(1)
  text(SR_dater[,"SSB"]/1000,SR_dater[,"recruit"]/1000,labels = rownames(SSB),adj = c(0,1),pos = 4,cex=0.8)  
mtext(side = 2,text = "Age-1 Recruitment (x1,000)",line = 4)



#Plot is looking like a beverton-holt would be our best bet - I also can't think of a reason overcompensation (cannibalism, density-dependent growth) is occuring.
###############
#Beverton-Holt#
###############
require("RTMB")
source(paste0(dat_path,"3_SR_func.R"))

dater<-list(SR_dater=SR_dater)
model_opt<-list(mod_type=2)
KMc_SR(params = list(alpha=1000,log_beta=-11,sigma=50))

###################
#Ricker - all data#
###################
dater<-list(SR_dater=SR_dater[as.character(seq(2009,2023,1)),])
model_opt<-list(mod_type=2)
obj<-MakeADFun(func = KMc_SR,parameters = list(alpha=1000,log_beta=-11,sigma=50))
opt<-nlminb(start = obj$par,objective =  obj$fn,gradient =  obj$gr)

sdr<-sdreport(obj = obj)
sdr

#Plot it out
par(mar=c(4.5,6,2,2))
plot(SR_dater[,"SSB"]/1000,SR_dater[,"recruit"]/1000,pch=21,col="black",bg="darkgrey",las=2,cex=1.2,
  xaxt="n",ylab="",xlab="SSB (x1,000 lbs)",xlim=c(0,600),ylim=c(0,1000))
  axis(1)
  text(SR_dater[,"SSB"]/1000,SR_dater[,"recruit"]/1000,labels = rownames(SSB),adj = c(0,1),pos = 4,cex=0.8)  
mtext(side = 2,text = "Age-1 Recruitment (x1,000)",line = 4)

pred_SSB<-seq(0,6000*1000,1000)
pred_recruit<-opt$par["alpha"]*pred_SSB*exp(-exp(opt$par["log_beta"])*pred_SSB)
lines(pred_SSB/1000,pred_recruit/1000,col="red")

#Then output the summary
#Fit a regression model to log(R/SSB) versus SSB and output the summary
Ricker.1b = lm(log(recruit/SSB)~SSB,data=SR_dater[as.character(seq(2008,2023,1)),])
summary(Ricker.1b)
#Extract the residual standard error for bias adjustment
RSE.Rick=sigma(Ricker.1b);

#Plot the unbiased predictions
Ricker.func<- function(SSB,alpha,beta){
  alpha*SSB*exp(-beta*SSB)
  } 
lines(pred_SSB/1000,
  Ricker.func(pred_SSB,exp(coef(Ricker.1b)[1]+0.5*RSE.Rick**2),abs(coef(Ricker.1b)[2]))/1000,lty=2,col="purple")


################
#B-H - all data#
################
dater<-list(SR_dater=SR_dater)
model_opt<-list(mod_type=1)
obj<-MakeADFun(func = KMc_SR,parameters = list(log_alpha=10,log_beta=10,sigma=50))
opt<-nlminb(start = obj$par,objective =  obj$fn,gradient =  obj$gr)

sdr<-sdreport(obj = obj)
sdr

#Plot it out
par(mar=c(4.5,6,2,2))
plot(SR_dater[,"SSB"]/1000,SR_dater[,"recruit"]/1000,pch=21,col="black",bg="darkgrey",las=2,cex=1.2,
  xaxt="n",ylab="",xlab="SSB (x1,000 lbs)",xlim=c(0,6000),ylim=c(0,1000))
  axis(1)
  text(SR_dater[,"SSB"]/1000,SR_dater[,"recruit"]/1000,labels = rownames(SSB),adj = c(0,1),pos = 4,cex=0.8)  
mtext(side = 2,text = "Age-1 Recruitment (x1,000)",line = 4)

pred_SSB<-seq(0,6000*1000,1000)
pred_recruit<-(exp(opt$par["log_alpha"])*pred_SSB)/(exp(opt$par["log_beta"])+pred_SSB)
lines(pred_SSB/1000,pred_recruit/1000,col="pink")

#Both Ricker and BH fire off into infinity - they appear to assume all the data is on assending limb

#####################
#BH with depensation#
#####################
dater<-list(SR_dater=SR_dater)
model_opt<-list(mod_type=3)
params_init<-list(log_alpha=5,log_beta=20,log_gamma=10,sigma=1)
obj<-MakeADFun(func = KMc_SR,parameters = params_init)

bounds_mat<-data.frame(param_name=names(x = params_init),lower=-Inf,upper=Inf)
  bounds_mat[bounds_mat[,"param_name"]=="log_gamma",c("lower","upper")]<-c(0,20)
opt_lower<-bounds_mat[match(x = names(obj$par), table=bounds_mat[,"param_name"]),"lower"]
opt_upper<-bounds_mat[match(x = names(obj$par), table=bounds_mat[,"param_name"]),"upper"]  

opt<-nlminb(start = obj$par,objective =  obj$fn,gradient =  obj$gr,control = list(eval.max=500,iter.max=1000),
  lower = opt_lower,upper = opt_upper)
opt
sdr<-sdreport(obj = obj)
sdr #Error on log(gamma) are kind of nonsensical

#Plot it out
par(mar=c(4.5,6,2,2))
plot(SR_dater[,"SSB"]/1000,SR_dater[,"recruit"]/1000,pch=21,col="black",bg="darkgrey",las=2,cex=1.2,
  xaxt="n",ylab="",xlab="SSB (x1,000 lbs)",xlim=c(0,600),ylim=c(0,1000))
  axis(1)
  text(SR_dater[,"SSB"]/1000,SR_dater[,"recruit"]/1000,labels = rownames(SSB),adj = c(0,1),pos = 4,cex=0.8)  
mtext(side = 2,text = "Age-1 Recruitment (x1,000)",line = 4)

pred_SSB<-seq(0,600000*1000,1000)
pred_recruit<-exp(opt$par["log_alpha"])*(pred_SSB^opt$par["gamma"])/(1+exp(opt$par["log_beta"])*(pred_SSB^opt$par["gamma"]))
lines(pred_SSB/1000,pred_recruit/1000,col="pink")

#########################
#Ricker with Depensation#
#########################
dater<-list(SR_dater=SR_dater)
model_opt<-list(mod_type=4)
params_init<-list(log_alpha=2,log_beta=-11,log_gamma=20,sigma=5)
obj<-MakeADFun(func = KMc_SR,parameters = params_init)


bounds_mat<-data.frame(param_name=names(x = params_init),lower=-Inf,upper=Inf)
  #bounds_mat[bounds_mat[,"param_name"]=="log_gamma",c("lower","upper")]<-c(5,Inf)
opt_lower<-bounds_mat[match(x = names(obj$par), table=bounds_mat[,"param_name"]),"lower"]
opt_upper<-bounds_mat[match(x = names(obj$par), table=bounds_mat[,"param_name"]),"upper"]  

opt<-nlminb(start = obj$par,objective =  obj$fn,gradient =  obj$gr,control = list(eval.max=500,iter.max=1000),
  lower = opt_lower, upper = opt_upper)
opt
sdr<-sdreport(obj = obj)
sdr #got it to fit...

#Plot it out
par(mar=c(4.5,6,2,2))
plot(SR_dater[,"SSB"]/1000,SR_dater[,"recruit"]/1000,pch=21,col="black",bg="darkgrey",las=2,cex=1.2,
  xaxt="n",ylab="",xlab="SSB (x1,000 lbs)",xlim=c(0,600),ylim=c(0,1000))
  axis(1)
  text(SR_dater[,"SSB"]/1000,SR_dater[,"recruit"]/1000,labels = rownames(SSB),adj = c(0,1),pos = 4,cex=0.8)  
mtext(side = 2,text = "Age-1 Recruitment (x1,000)",line = 4)

pred_SSB<-seq(0,6000*1000,1000)
(S/(S+exp(log_gamma)))*(exp(log_alpha)*S*exp(-exp(log_beta)*S))
opt$par["log_gamma"]

###########
#AR Ricker#
###########
source(paste0(dat_path,"3_SR_func.R"))

dater<-list(SR_dater=SR_dater)
KMc_SR_AR(params = list(log_alpha=10,log_beta=-11,log_gamma=2,sigma=50))

params_init<-list(log_alpha=10,log_beta=-11,log_gamma=2,sigma=50)
obj<-MakeADFun(func = KMc_SR_AR,parameters = params_init)
opt<-nlminb(start = obj$par,objective =  obj$fn,gradient =  obj$gr,control = list(eval.max=500,iter.max=1000))
opt
sdr<-sdreport(obj = obj)
sdr #Error on log(gamma) are kind of nonsensical

#Plot it out
par(mar=c(4.5,6,2,2))
plot(SR_dater[,"SSB"]/1000,SR_dater[,"recruit"]/1000,pch=21,col="black",bg="darkgrey",las=2,cex=1.2,
  xaxt="n",ylab="",xlab="SSB (x1,000 lbs)",xlim=c(0,600),ylim=c(0,1000))
  axis(1)
  text(SR_dater[,"SSB"]/1000,SR_dater[,"recruit"]/1000,labels = rownames(SSB),adj = c(0,1),pos = 4,cex=0.8)  
mtext(side = 2,text = "Age-1 Recruitment (x1,000)",line = 4)

pred_recruit<-numeric(length = nrow(SR_dater)-1)

for(i in 2:nrow(SR_dater)){
  pred_recruit[i-1]<-exp(opt$par["log_alpha"])*SR_dater[i,"SSB"]*exp(-exp(opt$par["log_beta"])*SR_dater[i,"SSB"]) + exp(opt$par["log_gamma"])*SR_dater[(i-1),"SSB"]
  }
  
#Pred_dater
pred_dater<-data.frame(SSB=SR_dater[2:nrow(SR_dater),"SSB"],pred_recruit=pred_recruit)
pred_dater<-pred_dater[order(pred_dater[,"SSB"]),]

lines(pred_dater[,"SSB"]/1000,pred_dater[,"pred_recruit"]/1000,col="red",type="b",pch=19)
#############################
#None of these really fit...#
#############################

#Okay going to try running a AR1 state space formulation of a ricker model...
require("rstan")

################
#AR formulation#
################
SRAR1_fit<-stan(file = paste0(dat_path,"3_SRAR1.stan"),model_name = "SRAR1",data = list(n=nrow(SR_dater),S=SR_dater[,"SSB"],log_R=log(SR_dater[,"recruit"])),
                iter=5000,chains=3,control = list(max_treedepth = 12,adapt_delta = 0.99))
SRAR1_fit
summary(SRAR1_fit)
traceplot(SRAR1_fit)

windows()
pairs(SRAR1_fit)#,pars = c("log_alpha","log_mean_R0"))


#########################
#Random walk formulation#
#########################
SRRW_fit<-stan(file = paste0(dat_path,"3_SRRW.stan"),model_name = "SRRW",data = list(n=nrow(SR_dater),S=SR_dater[,"SSB"],log_R=log(SR_dater[,"recruit"])),
                iter=20000,thin = 4,chains=3,control = list(max_treedepth = 13,adapt_delta = 0.99))
SRRW_fit
traceplot(SRRW_fit,pars = c("log_alpha","beta","phi_sd","eps_sd"))
traceplot(SRRW_fit,pars = c("phi"))
windows()
pairs(SRRW_fit,pars = c("log_alpha","beta","phi_sd","eps_sd","lp__"))

names(SRRW_fit)

plot(as.numeric(rownames(SR_dater)),
  exp(apply(X = as.matrix(SRRW_fit)[,paste0("log_alpha_w[",seq(1,25,1),"]")],MARGIN = 2,FUN = mean)),
  type="b")
















