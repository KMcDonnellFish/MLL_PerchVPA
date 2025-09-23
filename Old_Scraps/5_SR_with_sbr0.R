#Develop a random walk ricker SR model using the recK, R0, and sbr0 parameterization.

#Set Paths
dat_path<-"C:/Users/EU01237640/OneDrive - State of Minnesota - MN365/Mille_Lacs/YEP/PerchVPA/"

#I only have length/weight data from 1999 to present.  THey used scales prior to 1999 and I do not have those data.

#Load data
vpa_dat<-readRDS(file = paste0(dat_path,"opt_vpa.rds"))
LW_dat<-readRDS(file = paste0(dat_path,"LW_byage.rds"))
mat_dat<-readRDS(file = paste0(dat_path,"mat_byage.rds"))
sbr_dat<-readRDS(file = paste0(dat_path,"sbr_mat.rds"))

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

#Append the sbr0 data to this...
#make sure units of SSB and sbr0 are the same...
SR_dater<-cbind(SR_dater,sbr0=sbr_dat[rownames(sbr_dat)%in%rownames(SR_dater),"sbr0"]*2.20462) #Convert sbr0 from kg to lbs


plot(SR_dater[,"SSB"],SR_dater[,"recruit"],type="p",pch=21,col="darkgrey",bg="black")

#Need to provide some starting values
#Model
#Generate RW predictions

init_params<-expand.grid(log_recK=seq(0,100,0.5),log_R0=seq(0,10,0.5))
init_params$log_pred_R<-NA


for(j in 1:nrow(init_params)){
  n = nrow(SR_dater)
  S=SR_dater[,"SSB"]
  log_R=log(SR_dater[,"recruit"])
  sbr0=SR_dater[,"sbr0"]
  
   for(i in 1:n){
    if(i==1){
      log_pred_R<-rep(NA,nrow(SR_dater))
      w<-rep(NA,nrow(SR_dater))
      CR<-rep(NA,nrow(SR_dater))
      alpha<-rep(NA,nrow(SR_dater))
      beta<-rep(NA,nrow(SR_dater))
      
      w[1] = init_params[j,"log_recK"];
      
      } else {
        w[i] = w[(i-1)] + 0;
        }
    
    #Create easier intermediates
    CR[i] = exp(w[i]);
    alpha[i] = CR[i]/sbr0[i];
    beta[i] = log(CR[i])/(exp(init_params[j,"log_R0"])*sbr0[i]);
    #log scale Ricker    
    log_pred_R[i] = log(alpha[i]*S[i]*exp(-beta[i]*S[i]));
    }
  if(sum(is.infinite(log_pred_R))>0){
    init_params[j,"log_pred_R"]<-1
    } else {
      init_params[j,"log_pred_R"]<-0
      }
  }

init_params[init_params[,"log_pred_R"]==1,]
init_params[init_params[,"log_pred_R"]==0,]

colz<-c("grey","pink")

plot(init_params[,"log_R0"],init_params[,"log_recK"],type="p",pch=19,col=colz[init_params[,"log_pred_R"]+1])

#Initial values that don't provide Nans.
# log_recK = 100
# log_r0 = 6
  
mod_dater<-list(n=nrow(SR_dater),S=SR_dater[,"SSB"],log_R=log(SR_dater[,"recruit"]),sbr0=SR_dater[,"sbr0"])
SRRW_fit<-stan(file = paste0(dat_path,"5_SRRW_sbr0.stan"),model_name = "SRRW_sbr0",data = mod_dater,init = list(list(log_recK=50,log_R0=2)),
                iter=2500,thin = 1,chains=1,control = list(max_treedepth = 13,adapt_delta = 0.99),verbose = FALSE)

#Let's see how the CR/R0 formulation responds to different values
SR_sbro_fun<-function(S,sbr0,CR,R0){
  log_R<-log(S)-log(sbr0)+log(CR)*(1-(S/(R0*sbr0)))
  return(exp(log_R))
  }

curve(SR_sbro_fun(S = x,sbr0 = mean(SR_dater[,"sbr0"]),CR = 100,R0 = 10),from=0,to=50,ylab="Recruits",xlab="SSB",col="darkred",lwd=2,ylim=c(0,100))
  curve(SR_sbro_fun(S = x,sbr0 = mean(SR_dater[,"sbr0"]),CR = 50,R0 = 10),from=0,to=50,ylab="Recruits",xlab="SSB",col="red",lwd=2,add=TRUE)
  curve(SR_sbro_fun(S = x,sbr0 = mean(SR_dater[,"sbr0"]),CR = 25,R0 = 10),from=0,to=50,ylab="Recruits",xlab="SSB",col="pink",lwd=2,add=TRUE)

curve(SR_sbro_fun(S = x,sbr0 = mean(SR_dater[,"sbr0"]),CR = 100,R0 = 10),from=0,to=50,ylab="Recruits",xlab="SSB",col="darkred",lwd=2,ylim=c(0,100))
  curve(SR_sbro_fun(S = x,sbr0 = mean(SR_dater[,"sbr0"]),CR = 100,R0 = 15),from=0,to=50,ylab="Recruits",xlab="SSB",col="red",lwd=2,add=TRUE)
  curve(SR_sbro_fun(S = x,sbr0 = mean(SR_dater[,"sbr0"]),CR = 100,R0 = 5),from=0,to=50,ylab="Recruits",xlab="SSB",col="pink",lwd=2,add=TRUE)

  
  


mod_dater<-list(n=nrow(SR_dater),S=SR_dater[,"SSB"],log_R=log(SR_dater[,"recruit"]),sbr0=SR_dater[,"sbr0"])
SRRW_fit<-stan(file = paste0(dat_path,"5_SRRW_sbr0_V2.stan"),model_name = "SRRW_sbr0",data = mod_dater,
                iter=5000,thin = 1,chains=1)#,control = list(max_treedepth = 13,adapt_delta = 0.99),verbose = FALSE)

summary(SRRW_fit)

pairs(SRRW_fit,pars = c("log_recK","log_R0","eps_sd","phi_sd"))

traceplot(SRRW_fit,pars=c("log_recK","log_R0","eps_sd","phi_sd"))
traceplot(SRRW_fit,pars=c("log_recK","log_R0","eps_sd"))
traceplot(SRRW_fit,pars="phi")
stan_dens(object = SRRW_fit,pars = c("log_recK","log_R0","eps_sd","phi_sd"))
stan_dens(object = SRRW_fit,pars = c("log_recK","log_R0","eps_sd"))
plot(SRRW_fit,pars=c("log_recK","log_R0","eps_sd","phi_sd"))

#I cannot get any of these models to fit...  There's no compensation 
#Let's try fitting a depensation model in STAN.

mod_dater<-list(n=nrow(SR_dater),S=SR_dater[,"SSB"],log_R=log(SR_dater[,"recruit"]))
SR_ricker_depen_fit<-stan(file = paste0(dat_path,"5_SR_ricker_depen.stan"),model_name = "SR_ricker_depen",data = mod_dater,
                            iter=5000,thin = 1,chains=1)#,control = list(max_treedepth = 13,adapt_delta = 0.99),verbose = FALSE)
SR_ricker_depen_fit

summary(SR_ricker_depen_fit)
traceplot(SR_ricker_depen_fit)

#Beverton Holt depensation model
mod_dater<-list(n=nrow(SR_dater),S=SR_dater[,"SSB"],log_R=log(SR_dater[,"recruit"]))
SR_bh_depen_fit<-stan(file = paste0(dat_path,"5_SR_bh_depen.stan"),model_name = "SR_bh_depen",data = mod_dater,
                            iter=10000,thin = 1,chains=1)#,control = list(max_treedepth = 13,adapt_delta = 0.99),verbose = FALSE)
SR_bh_depen_fit



curve((x/(x+0))*x*exp(-0.68)*exp(-6.250861e-07*x),from=0,to=500000)
  curve((x/(x+50000))*x*exp(-0.68)*exp(-6.250861e-07*x),from=0,to=500000,col="darkred",add=TRUE)
  curve((x/(x+5000000))*x*exp(-0.68)*exp(-6.250861e-07*x),from=0,to=500000,col="red",add=TRUE)



windows()
  colz<-viridis::viridis(n = nrow(SR_dater))
  plot(SR_dater[,"SSB"],SR_dater[,"recruit"],pch=21,col="black",bg=colz)
    legend("topleft",legend=rownames(SR_dater),pch=21,col="black",pt.bg=colz)

curve(x/(x+250000),from=0,to=500000)
  


