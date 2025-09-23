#Run VPA code

######################
#Set Up data----
######################
#paths
base_path<-"C:/Users/EU01237640/OneDrive - State of Minnesota - MN365/Mille_Lacs/YEP/PerchVPA/"
dat_path<-paste0(base_path,"Data/")
code_pat<-paste0(base_path,"Scripts/")
out_path<-paste0(base_path,"Output/")

#Load data
#Catch/Survey data
catch_mat<-readRDS(paste0(dat_path,"catch_mat.rds"))
GN_mat<-readRDS(paste0(dat_path,"GN_mat.rds"))
offGN_mat<-readRDS(paste0(dat_path,"offGN_mat.rds"))
LW_array<-readRDS(paste0(dat_path,"LW_byage.rds"))
M<-readRDS(paste0(dat_path,"LorenzenM.rds"))

#VPA set up ----
ages <- as.numeric(gsub(pattern = "age",replacement = "",x = colnames(catch_mat)))
years <- as.numeric(rownames(catch_mat))


#Simple VPA----
#Based on code provided by MSU QFC.
#This code runs it for every terminal explotitation rate in seq(0,1,0.01)

ca <- catch_mat # catch at age for all gears

# NAs = 0
ca[which(is.na(ca))] <- 0

ny <- length(years)
na <- length(ages)
N <- matrix(0, nrow = ny, ncol = na) # numbers at age matrix

#Simulate for a bunch of different exploitation rates
Ut_valz<-seq(0,1,0.01)

N_array<-array(NA,dim = c(dim(N),length(Ut_valz)),dimnames = list(NULL,NULL,Ut_valz))

# STEP 1:
# initialize
for(i in 1:dim(N_array)[3]){
  Ut <- as.numeric(dimnames(N_array)[[3]][i])
  N[, na] <- ca[, na] / Ut
  N[!is.finite(N)] <- 0 # deal with divide by zero
  
  # STEP 2:
  # run backward VPA equation using Pope's approximation
  for (t in (ny - 1):1) {
    for (a in (na - 1):1) {
      N[t, a] <- N[t + 1, a + 1] * exp(M[a]) + ca[t, a] * exp(M[a] / 2)
    }
  }
  
  U_final <- ca / N
  # gill net vulnerability
  
  # STEP3:
  #reinitialize last year for U and N
  UbarT <- colMeans(U_final[1:(ny-1),]) # set for new Ubar for last year
  Vbar <- UbarT/UbarT[na-1] # vulnerability
  Vbar[na] <- 1
  N[ny,] <- ca[ny,] / (Vbar * Ut) # calculate N for last year
  U_final <- ca / N
  
  N_array[,,i]<-N
  
  }

#Quick plot
par(mar=c(4.5,6,2,2))
  plot(years,rowSums(N_array[,,"0.18"])/1000, type = "b", pch = 21, col = "black", bg = "darkgrey",ylim=c(0,50000),las=2,ylab="",xaxt="n")
  axis(1)
mtext(text = "Total Abundance (X1,000)",side = 2,line = 4)

#Most recent years
par(mar=c(4.5,6,2,2))
  plot(years,rowSums(N_array[,,"0.4"])/10000, type = "b", pch = 21, col = "black", bg = "darkgrey",ylim=c(0,200),las=2,ylab="",xaxt="n",xlim=c(2010,2024))
  axis(1)
mtext(text = "Total Abundance (X10,000)",side = 2,line = 4)


#Gulland Approach----

#This function takes the catch data, the terminal F, and age specific M and does a VPA.
Gulland_FUN<-function(catch,F_term,M,opt=FALSE){
  #Catch has to have form: Rows = Ages, Columns = Years
  #F_term is the terminal value of F
  #M is a vector of age specific values of M.
  #Opt is a logical to return the objective function value or the matrices
  #The objective function is the squared diff between Terminal F and F

  #Going to try my own implementation of the Gullands VPA with Pope's Approx.
  F_mat<-matrix(0,nrow = nrow(catch),ncol = ncol(catch),dimnames = dimnames(catch))
  N_mat<-matrix(0,nrow = nrow(catch),ncol = ncol(catch),dimnames = dimnames(catch))
  C_mat<-catch
  #Define Indices
  k<-nrow(F_mat) #ages
  l<-ncol(F_mat) #years
  
  #Inititalize using oldest age at last time step
  #1) Assume an Arbitrary F(k,l)
  F_mat[k,l]<-F_term
  #2) Calc. N[k,l]
  N_mat[k,l]<-((F_mat[k,l]+M[k])/F_mat[k,l])*C_mat[k,l]
  #3) Calc. N[k-1,l-1], N[k-2,l-2].... using Popes approx.
  for(i in 1:(k-1)){
    N_mat[k-i,l-i]<-N_mat[k-i+1,l-i+1]*exp(M[k-i])+C_mat[k-i,l-i]*exp(M[k-i]/2)
    F_mat[k-i,l-i]<- -log(N_mat[k-i+1,l-i+1]/N_mat[k-i,l-i]) - M[k-i]
    }
  #Use same process for all "complete" cohorts
  #4) Calc. F[k,l-1] by assuming F[k,l-1] = F[k-1,l-1]
  for(i in 1:(l-1)){
    F_mat[k,l-i]<-F_mat[k-1,l-i]
    N_mat[k,l-i]<-((F_mat[k,l-i]+M[k])/F_mat[k,l-i])*C_mat[k,l-i]
    for(j in 1:(k-1)){
      if( (k-j)>0 & (l-i-j)>0){
        N_mat[k-j,l-i-j]<-N_mat[k-j+1,l-i-j+1]*exp(M[k-j])+C_mat[k-j,l-i-j]*exp(M[k-j]/2)
        F_mat[k-j,l-i-j]<- -log(N_mat[k-j+1,l-i-j+1]/N_mat[k-j,l-i-j]) - M[k-j]
        }
      }
    }
  #Complete the "incomplete" cohorts
  for(j in 1:(k-1)){
    if((k-j)>0){
      F_mat[k-j,l]<-mean(F_mat[k-j,l-c(1,2,3)]) #Assume F[,l] is the mean of previous three years for each age
      N_mat[k-j,l]<-((F_mat[k-j,l]+M[k-j])/F_mat[k-j,l])*C_mat[k-j,l]*(1/(1-exp(-(F_mat[k-j,l]+M[k-j])))) #last year of VPN - assume catch is equal to N[a,l]*(F/Z)*(1-S)
      }
    for(i in 1:(k-1)){
      if((k-j-i)>0 & (l-i)>0){
        N_mat[k-j-i,l-i]<-N_mat[k-j-i+1,l-i+1]*exp(M[k-j-i])+C_mat[k-j-i,l-i]*exp(M[k-j-i]/2)
        F_mat[k-j-i,l-i]<- -log(N_mat[k-j-i+1,l-i+1]/N_mat[k-j-i,l-i]) - M[k-j-i]
        }
      }
    }
  
  if(opt==TRUE){
    return((F_mat[k,l]-F_mat[k-1,l])^2)
    } else {
      return(list(N_mat=N_mat,F_mat=F_mat))
      }
  }

opt_val<-nlminb(start = list(F_term = 0.5), #If initital value is <0.4, optimization will push terminal F to 0. 
        objective = Gulland_FUN, 
        catch = t(ca),
        M = M, 
        opt = TRUE,
        control = list(trace = 1, eval.max = 1000, iter.max = 1000),lower = 0.01,upper = Inf)
opt_val
opt_vpa<-Gulland_FUN(catch = t(ca),F_term = opt_val$par[1],M = M,opt = FALSE)

#Quick Look 
round(x = opt_vpa$N_mat,digits = 0)
round(opt_vpa$F_mat,digits = 4)

#Look at exploitation
round(t(ca)/opt_vpa$N_mat,digits = 4)
mean(round(t(ca)/opt_vpa$N_mat,digits = 4)[6,]) #0.354 is the average exploitation rate, use that for comparison plots


#Simple VPA vs Gulland Approach comparison plots----
#Total Abundance
plot(years,rowSums(N_array[,,"0.35"])/1000, type = "b", pch = 21, col = "steelblue", bg = "lightblue",ylim=c(0,35000),ylab="",las=2,xaxt="n")
  mtext(text = "Total Abund. (x1,000)",side = 2,line = 4.5)
  points(as.numeric(colnames(opt_vpa$N_mat)),colSums(opt_vpa$N_mat)/1000,type="b",pch=21,col="red",bg="pink")
  axis(1)
legend("topright",legend = c("Simple Mod","Gulland Approach"),pch=c(21,21),pt.bg=c("lightblue","pink"),col=c("steelblue","red"),lty=1)

plot(years,rowSums(N_array[,,"0.35"])/1000, type = "b", pch = 21, col = "steelblue", bg = "lightblue",ylim=c(0,5000),ylab="",las=2,xaxt="n",xlim=c(2010,2025))
  mtext(text = "Total Abund. (x1,000)",side = 2,line = 4.5)
  points(as.numeric(colnames(opt_vpa$N_mat)),colSums(opt_vpa$N_mat)/1000,type="b",pch=21,col="red",bg="pink")
  axis(1)
legend("topright",legend = c("Simple Mod","Gulland Approach"),pch=c(21,21),pt.bg=c("lightblue","pink"),col=c("steelblue","red"),lty=1)

#Age3+ Abundance
plot(years,rowSums(N_array[,3:6,"0.35"])/1000, type = "b", pch = 21, col = "steelblue", bg = "lightblue",ylim=c(0,8000),ylab="Abund. 3+ (x1,000)",main="VPA Est. Age3+ YEP",las=2,xaxt="n")
  axis(1)
  points(as.numeric(colnames(opt_vpa$N_mat)),colSums(opt_vpa$N_mat[3:6,])/1000,type="b",pch=21,col="red",bg="pink")
legend("topright",legend = c("Simple Mod","Gulland Approach"),pch=c(21,21),pt.bg=c("lightblue","pink"),col=c("steelblue","red"),lty=1)

plot(years,rowSums(N_array[,3:6,"0.35"])/1000, type = "b", pch = 21, col = "steelblue", bg = "lightblue",ylim=c(0,2000),ylab="Abund. 3+ (x1,000)",main="VPA Est. Age3+ YEP",las=2,xaxt="n",xlim=c(2010,2025))
  axis(1)
  points(as.numeric(colnames(opt_vpa$N_mat)),colSums(opt_vpa$N_mat[3:6,])/1000,type="b",pch=21,col="red",bg="pink")
legend("topright",legend = c("Simple Mod","Gulland Approach"),pch=c(21,21),pt.bg=c("lightblue","pink"),col=c("steelblue","red"),lty=1)

#Age 1 Abundance
plot(years,N_array[,1,"0.35"]/1000, type = "b", pch = 21, col = "steelblue", bg = "lightblue",ylim=c(0,30000),ylab="",main="VPA Est. Age1 YEP",las=2,xaxt="n")
  points(as.numeric(colnames(opt_vpa$N_mat)),opt_vpa$N_mat[1,]/1000,type="b",pch=21,col="red",bg="pink")
  axis(1)
  mtext(text = "Abund. Age1 (x1,000)",side = 2,line = 4.5)
legend("topright",legend = c("Simple Mod","Gulland Approach"),pch=c(21,21),pt.bg=c("lightblue","pink"),col=c("steelblue","red"),lty=1)

plot(years,N_array[,1,"0.35"]/1000, type = "b", pch = 21, col = "steelblue", bg = "lightblue",ylim=c(0,2500),ylab="",main="VPA Est. Age1 YEP",las=2,xaxt="n",xlim=c(2010,2025))
  points(as.numeric(colnames(opt_vpa$N_mat)),opt_vpa$N_mat[1,]/1000,type="b",pch=21,col="red",bg="pink")
  axis(1)
  mtext(text = "Abund. Age1 (x1,000)",side = 2,line = 4.5)
legend("topright",legend = c("Simple Mod","Gulland Approach"),pch=c(21,21),pt.bg=c("lightblue","pink"),col=c("steelblue","red"),lty=1)

#Plot Age5/6 F
# plot(NA,NA, type = "n", pch = 21, col = "red", bg = "pink",xlim=c(min(years),max(years)),ylim=c(0,2),ylab="Z",main="Age 5 and 6",las=2,xaxt="n")
#   axis(1)
#   polygon(c(years,rev(years)),c(opt_vpa$F_mat[6,]+M,rep(0,length(years))),col = "pink",border=NA)
#   polygon(c(years,rev(years)),c(rep(M,length(years)),rep(0,length(years))),col = "darkgrey",border = NA)
# legend("topright",legend = c("F","M"),fill=c("pink","darkgrey"))


#Mortality Plots----
#Plot as a barplot instaed
par(mfcol=c(3,2))
for(i in 1:6){
  x<-barplot(rbind(M[i],opt_vpa$F_mat[i,]),col=c("darkgrey","pink"),xaxt="n",las=2,ylim=c(0,2),ylab="Z")
    axis(1,at = x[years%in%seq(1990,2020,5)],labels = years[years%in%seq(1990,2020,5)])
    legend("topright",legend = c("F","M"),fill=c("pink","darkgrey"),title = paste0("Age ",i))
  }

#Calculate F for simple VPA - double check this....
#C/N = U = F/(F+M)
#F = UM/(1-U)
simp_Ut<-t(ca/N_array[,,"0.35"])
M_mat<-matrix(M,nrow = dim(simp_Ut)[1],ncol = dim(simp_Ut)[2])
simp_F<-(simp_Ut*M_mat)/(1-simp_Ut)

par(mfcol=c(3,2))
for(i in 1:6){
  x<-barplot(rbind(M[i],simp_F[i,]),col=c("darkgrey","pink"),xaxt="n",las=2,ylim=c(0,2),ylab="Z")
    axis(1,at = x[years%in%seq(1990,2020,5)],labels = years[years%in%seq(1990,2020,5)])
    legend("topright",legend = c("F","M"),fill=c("pink","darkgrey"),title = paste0("Age ",i))
  }

#Exploitation plots----
#Exploitation rate by age
colz<-c("#fde725","#5ec962","#21918c","#3b528b","#d093e0ff","#440154")
gulland_Ut<-t(ca)/opt_vpa$N_mat
plot(0,0,type = "n",xlim = summary(as.numeric(colnames(gulland_Ut)))[c("Min.","Max.")],ylim=c(0,1),las=2,xaxt="n",ylab="Exploitation Rate",xlab="Year")
  axis(1)
  for(i in 1:nrow(gulland_Ut)){
    lines(as.numeric(colnames(gulland_Ut)), gulland_Ut[i,],lwd=2,col=colz[i])
    points(as.numeric(colnames(gulland_Ut)), gulland_Ut[i,],pch=21,col="darkgrey",bg=colz[i])
    }
legend("topright",legend=paste0("age ",seq(1,6,1)),pch=21,col="darkgrey",pt.bg=colz,pt.cex=2)

colz<-c("#fde725","#5ec962","#21918c","#3b528b","#d093e0ff","#440154")
par(mfcol=c(1,1))
plot(0,0,type = "n",xlim = summary(as.numeric(colnames(simp_Ut)))[c("Min.","Max.")],ylim=c(0,1),las=2,xaxt="n",ylab="Exploitation Rate",xlab="Year")
  axis(1)
  for(i in 1:nrow(simp_Ut)){
    lines(as.numeric(colnames(simp_Ut)), simp_Ut[i,],lwd=2,col=colz[i])
    points(as.numeric(colnames(simp_Ut)), simp_Ut[i,],pch=21,col="darkgrey",bg=colz[i])
    }
legend("topright",legend=paste0("age ",seq(1,6,1)),pch=21,col="darkgrey",pt.bg=colz,pt.cex=2)


#Abundance estimates vs survey CPUE----
#How does the age3+ plot compare to the nets CPUE?
#Age3+ Abundance
par(mar=c(4.5,6,2,6))
plot(years,rowSums(N_array[,3:6,"0.35"])/1000, type = "b", pch = 21, col = "steelblue", bg = "lightblue",ylim=c(0,7000),ylab="Abund. 3+ (x1,000)",main="VPA Est. Age3+ YEP",las=2,xaxt="n")
  points(as.numeric(colnames(opt_vpa$N_mat)),colSums(opt_vpa$N_mat[3:6,])/1000,type="b",pch=21,col="red",bg="pink")
  axis(1)
  
  par(new=TRUE)
  plot(years,rowSums(GN_mat[-nrow(GN_mat),3:6]), type = "b", pch = 22, col = "black", bg = "darkgrey",ylim=c(0,150),ylab="",main="",yaxt="n",xaxt="n")
  mtext("Nearshore GN CPUE (#Age3+/Net)",side = 4,line = 3.5)
  axis(4,las=2)
legend("topright",legend = c("Simple Mod","Gulland Approach","GN CPUE"),pch=c(21,21,22),pt.bg=c("lightblue","pink","darkgrey"),col=c("steelblue","red","black"),lty=1)


par(mar=c(4.5,6,2,6))
plot(years,rowSums(N_array[,3:6,"0.42"])/1000, type = "b", pch = 21, col = "steelblue", bg = "lightblue",ylim=c(0,7000),ylab="Abund. 3+ (x1,000)",main="VPA Est. Age3+ YEP",las=2,xaxt="n")
  points(as.numeric(colnames(opt_vpa$N_mat)),colSums(opt_vpa$N_mat[3:6,])/1000,type="b",pch=21,col="red",bg="pink")
  axis(1)
  
  par(new=TRUE)
  plot(years[years>=1999],rowSums(offGN_mat[-nrow(offGN_mat),3:6]), type = "b", pch = 22, col = "black", bg = "darkgrey",xlim=c(min(years),max(years)),ylim=c(0,100),ylab="",main="",yaxt="n",xaxt="n")
  mtext("Offshore GN CPUE (#Age3+/Net)",side = 4,line = 3.5)
  axis(4,las=2)
legend("topright",legend = c("Simple Mod","Gulland Approach","GN CPUE"),pch=c(21,21,22),pt.bg=c("lightblue","pink","darkgrey"),col=c("steelblue","red","black"),lty=1)

#Exploitation of SSB----
#Let's look at Ut of SSB
#Pull in weight at age data and maturity data
LW_dat<-readRDS(file = paste0(dat_path,"LW_byage.rds"))
mat_dat<-readRDS(file = paste0(dat_path,"mat_byage.rds"))

#SSB through time - only have these data from 1999 to present for now.  Older data is squirreled away in a ancient database apparently.

#I'm going backcast weight-at-age and maturity at age from 1999 to 1985.  Assume weight-at-age, and maturity-at-age from 1986-1998 = Mean values from 1999-2004
hist_weights<-apply(LW_dat[as.character(1999:2004),,"weight"],MARGIN = 2,FUN = mean)
hist_mat<-apply(mat_dat[as.character(1999:2004),],MARGIN = 2,FUN = mean)

hist_weights<-matrix(hist_weights,nrow = length(1986:1998),ncol=dim(LW_dat)[2],byrow=TRUE,dimnames=list(seq(1986,1998),dimnames(LW_dat)[[2]]))
hist_mat<-matrix(hist_mat,nrow = length(1986:1998),ncol=dim(mat_dat)[2],byrow=TRUE,dimnames=list(seq(1986,1998),dimnames(mat_dat)[[2]]))

all_weights<-rbind(hist_weights,LW_dat[,,"weight"])
all_mat<-rbind(hist_mat,mat_dat)

SSB_mat<-opt_vpa$N_mat*t(all_weights[,as.character(1:6)])*0.00220462*t(all_mat[,as.character(1:6)]) #1g = 0.00220462lbs

plot(as.numeric(colnames(SSB_mat)),apply(SSB_mat,MARGIN = 2,FUN = sum)/1000,type="b",lwd=1.5,pch=21,col="black",bg="darkgrey",las=2,xaxt="n",
      xlab="Year",ylab="SSB (x1,000 lbs)")
axis(1)

#Exploitation of SSB
catch_weight<-t(catch_mat*all_weights[,as.character(1:6)]*0.00220462)
U_ssb<-apply(X = catch_weight,MARGIN = 2,FUN = sum)/apply(X = SSB_mat,MARGIN = 2,FUN = sum)

plot(as.numeric(names(U_ssb)),U_ssb,type="b",lwd=1.5,pch=21,col="black",bg="darkgrey",las=2,xaxt="n",ylim=c(0,0.4),
      xlab="Year",ylab="Ut of SSB (Yield/SSB)")
axis(1)

#Plot SSB and U together.
par(mar=c(4.5,5,2,5))
  plot(as.numeric(colnames(SSB_mat)),apply(SSB_mat,MARGIN = 2,FUN = sum)/1000,type="b",lwd=1.5,pch=21,col="black",bg="darkgrey",las=2,xaxt="n",
      xlab="Year",ylab="")
  mtext(text = "SSB (x1,000 lbs)",side = 2,line = 3.5)
  axis(1)
  par(new=TRUE)
  plot(as.numeric(names(U_ssb)),U_ssb,type="b",lwd=1.5,pch=22,col="red",bg="pink",las=2,xaxt="n",
      xlab="Year",ylab="",yaxt="n")
  axis(4,las=2)
  mtext(text = "Ut of SSB (Yield/SSB)",side = 4,line = 3.5)
legend("topright",legend=c("SSB","Ut of SSB"),col=c("black","red"),pch=c(21,22),pt.bg=c("darkgrey","pink"))

#Scatter Plot of SSB vs Ut of SSB
par(mar=c(4.5,4,2,2))
plot(apply(SSB_mat,MARGIN = 2,FUN = sum)/1000,U_ssb,pch=21,col="black",bg="darkgrey",cex=1.5,xlab="SSB (x1,000 lbs)",ylab="Ut of SSB (Yield/SSB)",las=2,xaxt="n")
axis(1)

#Plot SSB and total U together
par(mar=c(4.5,5,2,5))
  plot(as.numeric(colnames(SSB_mat)),apply(SSB_mat,MARGIN = 2,FUN = sum)/1000,type="b",lwd=1.5,pch=21,col="black",bg="darkgrey",las=2,xaxt="n",
      xlab="Year",ylab="")
  mtext(text = "SSB (x1,000 lbs)",side = 2,line = 3.5)
  axis(1)
  par(new=TRUE)

  Ut<-t(apply(X = catch_mat,MARGIN = 1,FUN = sum))/apply(X = opt_vpa$N_mat,MARGIN = 2,FUN = sum)

  plot(as.numeric(colnames(Ut)),Ut,type="b",lwd=1.5,pch=22,col="red",bg="pink",las=2,xaxt="n",
      xlab="Year",ylab="",yaxt="n")
  axis(4,las=2)
  mtext(text = "Ut (Catch/N)",side = 4,line = 3.5)
legend("topright",legend=c("SSB","Ut"),col=c("black","red"),pch=c(21,22),pt.bg=c("darkgrey","pink"))

#Scatter Plot of SSB vs Ut of SSB
par(mar=c(4.5,4,2,2))
plot(apply(SSB_mat,MARGIN = 2,FUN = sum)/1000,Ut,pch=21,col="black",bg="darkgrey",cex=1.5,xlab="SSB (x1,000 lbs)",ylab="Ut (Catch/N)",las=2,xaxt="n")
axis(1)


#Gulland VPA sensitivity----
#######################################
#Sensitivity plots for the Gulland Mod#
#######################################
#How does Age3+ abundance move with different values of Fterm?

##Terminal F----

F_termz_abund<-matrix(NA,nrow=6,ncol = length(years))
rownames(F_termz_abund)<-c(0.05,0.25,opt_val$par,0.75,1,1.5)
colnames(F_termz_abund)<-years
F_termz_Fval<-F_termz_abund


for(z in 1:nrow(F_termz_abund)){
  vpa_valz<-Gulland_FUN(catch = t(ca),F_term = as.numeric(rownames(F_termz_abund))[z],M = M,opt = FALSE)  
  F_termz_abund[z,]<-colSums(vpa_valz$N_mat[3:6,])
  F_termz_Fval[z,]<-vpa_valz$F_mat[6,]
  }

#Plot abundance
colz<-c("#f72585", "#b5179e", "#7209b7", "#3a0ca3", "#4361ee", "#4cc9f0")
plot(years,F_termz_abund[1,]/1000,type="l",las=2,xaxt="n",ylab="Age 3+ Abund. (x1,000)",lwd=1,col=colz[1])
  for(i in 2:nrow(F_termz_abund)){
    lines(years,F_termz_abund[i,]/1000,lwd=1,col=colz[i])
    }
  axis(1)
legend("topright",legend=round(as.numeric(rownames(F_termz_abund)),digits = 3),col=colz,lwd=2,title = "Values of F_term")
  
colz<-c("#f72585", "#b5179e", "#7209b7", "#3a0ca3", "#4361ee", "#4cc9f0")
plot(years,F_termz_abund[1,]/1000,type="l",las=2,xaxt="n",ylab="Age 3+ Abund. (x1,000)",lwd=1,col=colz[1],xlim=c(2010,2025),ylim=c(0,500))
  for(i in 2:nrow(F_termz_abund)){
    lines(years,F_termz_abund[i,]/1000,lwd=1,col=colz[i])
    }
  axis(1)
legend("topright",legend=round(as.numeric(rownames(F_termz_abund)),digits = 3),col=colz,lwd=2,title = "Values of F_term")

#Plot Age 6 F
#Plot abundance
colz<-c("#f72585", "#b5179e", "#7209b7", "#3a0ca3", "#4361ee", "#4cc9f0")
par(mar=c(4.5,6,2,2))
plot(years,F_termz_Fval[1,],type="l",las=2,xaxt="n",ylab="",lwd=1,col=colz[1],ylim=c(0,2.5))
  for(i in 2:nrow(F_termz_Fval)){
    lines(years,F_termz_Fval[i,],lwd=1,col=colz[i])
    }
  axis(1)
  mtext(side = 2,text = "Age6 F",line = 4.5)
legend("topright",legend=round(as.numeric(rownames(F_termz_Fval)),digits = 3),col=colz,lwd=2,title = "Values of F_term")

colz<-c("#f72585", "#b5179e", "#7209b7", "#3a0ca3", "#4361ee", "#4cc9f0")  
par(mar=c(4.5,6,2,2))
plot(years,F_termz_Fval[1,],type="l",las=2,xaxt="n",ylab="",lwd=1,col=colz[1],ylim=c(0,2.5),xlim=c(2010,2025))
  for(i in 2:nrow(F_termz_Fval)){
    lines(years,F_termz_Fval[i,],lwd=1,col=colz[i])
    }
  axis(1)
  mtext(side = 2,text = "Age6 F",line = 4.5)
legend("topright",legend=round(as.numeric(rownames(F_termz_Fval)),digits = 3),col=colz,lwd=2,title = "Values of F_term")

#VPA is not sensitive to value of F_term as long as F_term doesn't approach 0.
#Values of F seem to converge in after about 3 years regardless of starting point

##############
#5 year Retro#
##############
##Retrospective----

retro_abund<-matrix(NA,nrow=6,ncol = length(years))
rownames(retro_abund)<-seq(0,5)
colnames(retro_abund)<-years
retro_F<-retro_abund

for(z in 1:nrow(retro_abund)){
  sim_years<-years[1:(length(years)-as.numeric(rownames(retro_abund)[z]))]
  opt_val<-nlminb(start = list(F_term = 2), #If initital value is <0.4, optimization will push terminal F to 0. 
        objective = Gulland_FUN, 
        catch = t(ca)[,as.character(sim_years)],
        M = M, 
        opt = TRUE,
        control = list(trace = 1, eval.max = 1000, iter.max = 1000),lower = 0.1,upper = Inf)
  opt_vpa<-Gulland_FUN(catch = t(ca)[,as.character(sim_years)],F_term = opt_val$par,M = M,opt = FALSE)
  retro_abund[z,colnames(retro_abund)%in%sim_years]<-colSums(opt_vpa$N_mat[3:6,])
  retro_F[z,colnames(retro_F)%in%sim_years]<-opt_vpa$F_mat[6,]
  }

#Plot abundance
colz<-c("#f72585", "#b5179e", "#7209b7", "#3a0ca3", "#4361ee", "#4cc9f0")
plot(years,retro_abund[1,]/1000,type="l",las=2,xaxt="n",ylab="Age 3+ Abund. (x1,000)",lwd=2,col=colz[1])
  for(i in 2:nrow(retro_abund)){
    lines(years,retro_abund[i,]/1000,lwd=2,col=colz[i])
    }
  axis(1)
legend("topright",legend=rownames(retro_abund),col=colz,lwd=2,title = "Years removed from analysis")

colz<-c("#f72585", "#b5179e", "#7209b7", "#3a0ca3", "#4361ee", "#4cc9f0")
plot(years,retro_abund[1,]/1000,type="l",las=2,xaxt="n",ylab="Age 3+ Abund. (x1,000)",lwd=2,col=colz[1],xlim=c(2010,2025),ylim=c(0,500))
  for(i in 2:nrow(retro_abund)){
    lines(years,retro_abund[i,]/1000,lwd=2,col=colz[i])
    }
  axis(1)
legend("topright",legend=rownames(retro_abund),col=colz,lwd=2,title = "Years removed from analysis")

#Plot Age6 F values
colz<-c("#f72585", "#b5179e", "#7209b7", "#3a0ca3", "#4361ee", "#4cc9f0")
plot(years,retro_F[1,],type="l",las=2,xaxt="n",ylab="Age6 F",lwd=2,col=colz[1])
  for(i in 2:nrow(retro_F)){
    lines(years,retro_F[i,],lwd=2,col=colz[i])
    }
  axis(1)
legend("topright",legend=rownames(retro_F),col=colz,lwd=2,title = "Years removed from analysis")

colz<-c("#f72585", "#b5179e", "#7209b7", "#3a0ca3", "#4361ee", "#4cc9f0")
plot(years,retro_F[1,],type="l",las=2,xaxt="n",ylab="Age6 F",lwd=2,col=colz[1],xlim=c(2010,2025))
  for(i in 2:nrow(retro_F)){
    lines(years,retro_F[i,],lwd=2,col=colz[i])
    }
  axis(1)
legend("topright",legend=rownames(retro_F),col=colz,lwd=2,title = "Years removed from analysis")

#Save final output----
#Write the VPA abundance as a .rds
opt_val<-nlminb(start = list(F_term = 2), #If initital value is <0.4, optimization will push terminal F to 0. 
        objective = Gulland_FUN, 
        catch = t(ca), 
        opt = TRUE,
        M = M,
        control = list(trace = 1, eval.max = 1000, iter.max = 1000),lower = 0.1,upper = Inf)
opt_vpa<-Gulland_FUN(catch = t(ca),F_term = opt_val$par,M = M,opt = FALSE)
saveRDS(object = opt_vpa,file = paste0(dat_path,"opt_vpa.rds"))