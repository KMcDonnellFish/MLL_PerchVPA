######################
#Set Up data----
######################

#paths
dat_path<-"C:/Users/EU01237640/OneDrive - State of Minnesota - MN365/Mille_Lacs/YEP/PerchVPA/"

#Pull data out of the .dat file
datfile<-readLines(con = paste0(dat_path,"Perch2025.dat"))

#Grab VPA yearz
first_year<-as.numeric(datfile[which(grepl(pattern = "# First year of harvest data",x = datfile))+1])
last_year<-as.numeric(datfile[which(grepl(pattern = "# Last year of harvest data used in the analysis",x = datfile))+1])

#Pull in Harvest by year
harv_index_start<-which(grepl(pattern = "# Mille Lacs Yellow Perch Harvest Data",x = datfile))+2
harv_index_end<-harv_index_start+(last_year-first_year)

catch_dat<-datfile[harv_index_start:harv_index_end]
for(i in 1:length(catch_dat)){
  if(i==1){
    catch_mat<-matrix(NA,nrow = length(catch_dat),ncol = 6,dimnames = list(seq(first_year,last_year),paste0("age",seq(1,6))))
    }
  if(grepl(pattern = "\t",x = catch_dat[i])) catch_mat[i,]<-as.numeric(unlist(strsplit(x =catch_dat[i],split = "\t")))
  if(grepl(pattern = " ",x = catch_dat[i])){
    x<-as.numeric(unlist(strsplit(x =catch_dat[i],split = " ")))
    catch_mat[i,]<-x[!is.na(x)]
    }
  }
#Can ignore the NAs

#Get net CPUE data. from dat file
GN_index_start<-which(grepl(pattern = "# Mille Lacs Standard Gill Net Data",x = datfile))+2
GN_index_end<-GN_index_start+(last_year-first_year)+1 #last year is the 2024 GN season that matches with the unavailable 2025 fishing season
GN_dat<-datfile[GN_index_start:GN_index_end]
for(i in 1:length(GN_dat)){
  if(i==1){
    GN_mat<-matrix(NA,nrow = length(GN_dat),ncol = 6,dimnames = list(seq(first_year-1,last_year),paste0("age",seq(1,6))))
    }
  if(grepl(pattern = "\t",x = GN_dat[i])) GN_mat[i,]<-as.numeric(unlist(strsplit(x =GN_dat[i],split = "\t")))
  if(grepl(pattern = " ",x = GN_dat[i])){
    x<-as.numeric(unlist(strsplit(x = GN_dat[i],split = " ")))
    GN_mat[i,]<-x[!is.na(x)]
    }
  #print(i)
  }

#Get net CPUE data. from dat file
GN_index_start<-which(grepl(pattern = "# Mille Lacs Offshore Gill Net Data",x = datfile))+1
GN_index_end<-GN_index_start+(last_year-1998) #last year is the 2024 GN season that matches with the unavailable 2025 fishing season
offGN_dat<-datfile[GN_index_start:GN_index_end]

for(i in 1:length(offGN_dat)){
  if(i==1){
    offGN_mat<-matrix(NA,nrow = length(offGN_dat),ncol = 6,dimnames = list(seq(1998,last_year),paste0("age",seq(1,6))))
    }
  if(grepl(pattern = "\t",x = offGN_dat[i])) offGN_mat[i,]<-as.numeric(unlist(strsplit(x =offGN_dat[i],split = "\t")))
  if(grepl(pattern = " ",x = offGN_dat[i])){
    x<-as.numeric(unlist(strsplit(x = offGN_dat[i],split = " ")))
    offGN_mat[i,]<-x[!is.na(x)]
    }
  #print(i)
  }

##Summary Plots of input data ----
#find colors for ages
colz<-c("#fde725","#5ec962","#21918c","#3b528b","#d093e0ff","#440154")

catch_plot<-apply(X = catch_mat,MARGIN = 1,FUN = cumsum)
windows()
  par(mar=c(4.5,7,2,2))
  plot(0,0,type="n",xlim=c(1985,2025),ylim=c(0,1200000),cex.axis=1.2,las=2,xaxt="n",xlab="Year",ylab="",cex.lab=1.2)
  axis(1,cex.axis=1.2)
  plot_yearz<-as.numeric(dimnames(catch_plot)[[2]])
  for(i in nrow(catch_plot):1){
    polygon(x = c(plot_yearz,rev(plot_yearz)),y = c(catch_plot[i,],rep(0,length(catch_plot[i,]))),col=colz[i],border=NA)
    }
  mtext(text = "# YEP Harvested",side = 2,line = 5.5,cex=1.2)
legend("topright",legend=paste0("age",seq(1,6)),fill=colz,cex=1.2)

GN_plot<-apply(X = GN_mat,MARGIN = 1,FUN = cumsum)
windows()
  par(mar=c(4.5,7,2,2))
  plot(0,0,type="n",xlim=c(1985,2025),ylim=c(0,250),cex.axis=1.2,las=2,xaxt="n",xlab="Year",ylab="",cex.lab=1.2)
  axis(1,cex.axis=1.2)
  plot_yearz<-as.numeric(dimnames(GN_plot)[[2]])
  for(i in nrow(GN_plot):1){
    polygon(x = c(plot_yearz,rev(plot_yearz)),y = c(GN_plot[i,],rep(0,length(GN_plot[i,]))),col=colz[i],border=NA)
    }
  mtext(text = "CPUE Nearshore Fall Survey",side = 2,line = 5.5,cex=1.2)
legend("topright",legend=paste0("age",seq(1,6)),fill=colz,cex=1.2)



#Lets do a polygon plot of catch at age
colnames(offGN_mat)

 



























#VPA set up
## Setup ####
ages <- as.numeric(gsub(pattern = "age",replacement = "",x = colnames(catch_mat)))
years <- as.numeric(rownames(catch_mat))
#wbar <- colMeans(dat$wa)
ca <- catch_mat # catch at age for all gears

# NAs = 0
ca[which(is.na(ca))] <- 0

# Inputs
M <- 0.2
S <- exp(-M)

ny <- length(years)
na <- length(ages)
N <- matrix(0, nrow = ny, ncol = na) # numbers at age matrix

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
      N[t, a] <- N[t + 1, a + 1] * exp(M) + ca[t, a] * exp(M / 2)
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

#Quick plots
plot(as.numeric(dimnames(N_array)[[3]]),colSums(N_array[dim(N_array)[1],,]))
plot(years,rowSums(N_array[,,"0.42"]), type = "b", pch = 21, col = "black", bg = "darkgrey",ylim=c(0,5500000))

Gulland_FUN<-function(catch,F_term,opt=FALSE){
  #Catch has to have form: Rows = Ages, Columns = Years
  #F_term is the terminal value of F
  #Opt is a logical to return the objective function value or the matrices
  #The objective function is the squared diff between Terminal F and F
  
  #Going to try my own implementation of the Gullands VPA with Pope's Approx.
  F_mat<-matrix(0,nrow = nrow(catch),ncol = ncol(catch),dimnames = dimnames(catch))
  N_mat<-matrix(0,nrow = nrow(catch),ncol = ncol(catch),dimnames = dimnames(catch))
  C_mat<-catch
  #Define Indices
  k<-nrow(F_mat)
  l<-ncol(F_mat)
  
  #Inititalize using oldest age at last time step
  #1) Assume an Arbitrary F(k,l)
  F_mat[k,l]<-F_term
  #2) Calc. N[k,l]
  N_mat[k,l]<-((F_mat[k,l]+M)/F_mat[k,l])*C_mat[k,l]
  #3) Calc. N[k-1,l-1], N[k-2,l-2].... using Popes approx.
  for(i in 1:(k-1)){
    N_mat[k-i,l-i]<-N_mat[k-i+1,l-i+1]*exp(M)+C_mat[k-i,l-i]*exp(M/2)
    F_mat[k-i,l-i]<- -log(N_mat[k-i+1,l-i+1]/N_mat[k-i,l-i]) - M
    }
  #Use same process for all "complete" cohorts
  #4) Calc. F[k,l-1] by assuming F[k,l-1] = F[k-1,l-1]
  for(i in 1:(l-1)){
    F_mat[k,l-i]<-F_mat[k-1,l-i]
    N_mat[k,l-i]<-((F_mat[k,l-i]+M)/F_mat[k,l-i])*C_mat[k,l-i]
    for(j in 1:(k-1)){
      if( (k-j)>0 & (l-i-j)>0){
        N_mat[k-j,l-i-j]<-N_mat[k-j+1,l-i-j+1]*exp(M)+C_mat[k-j,l-i-j]*exp(M/2)
        F_mat[k-j,l-i-j]<- -log(N_mat[k-j+1,l-i-j+1]/N_mat[k-j,l-i-j]) - M
        }
      }
    }
  #Complete the "incomplete" cohorts
  for(j in 1:(k-1)){
    if((k-j)>0){
      F_mat[k-j,l]<-mean(F_mat[k-j,l-c(1,2,3)]) #Assume F[,l] is the mean of previous three years for each age
      N_mat[k-j,l]<-((F_mat[k-j,l]+M)/F_mat[k-j,l])*C_mat[k-j,l]*(1/(1-exp(-(F_mat[k-j,l]+M)))) #last year of VPN - assume catch is equal to N[a,l]*(F/Z)*(1-S)
      }
    for(i in 1:(k-1)){
      if((k-j-i)>0 & (l-i)>0){
        N_mat[k-j-i,l-i]<-N_mat[k-j-i+1,l-i+1]*exp(M)+C_mat[k-j-i,l-i]*exp(M/2)
        F_mat[k-j-i,l-i]<- -log(N_mat[k-j-i+1,l-i+1]/N_mat[k-j-i,l-i]) - M
        }
      }
    }
  
  if(opt==TRUE){
    return((F_mat[k,l]-F_mat[k-1,l])^2)
    } else {
      return(list(N_mat=N_mat,F_mat=F_mat))
      }
  }

opt_val<-nlminb(start = list(F_term = 2), #If initital value is <0.4, optimization will push terminal F to 0. 
        objective = Gulland_FUN, 
        catch = t(ca), 
        opt = TRUE,
        control = list(trace = 1, eval.max = 1000, iter.max = 1000),lower = 0.1,upper = Inf)
opt_vpa<-Gulland_FUN(catch = t(ca),F_term = opt_val$par,opt = FALSE)

#Total Abundance
plot(years,rowSums(N_array[,,"0.42"])/1000, type = "b", pch = 21, col = "steelblue", bg = "lightblue",ylim=c(0,5500),ylab="Total Abund. (x1,000)",las=2,xaxt="n")
  points(as.numeric(colnames(opt_vpa$N_mat)),colSums(opt_vpa$N_mat)/1000,type="b",pch=21,col="red",bg="pink")
  axis(1)
legend("topright",legend = c("Simple Mod","Gulland Approach"),pch=c(21,21),pt.bg=c("lightblue","pink"),col=c("steelblue","red"),lty=1)

#Total Abundance - plot 2010 - present
plot(years[years>=2010],rowSums(N_array[years>=2010,,"0.42"])/1000, type = "b", pch = 21, col = "steelblue", bg = "lightblue",ylim=c(0,1000),ylab="Total Abund. (x1,000)",las=2,xaxt="n")
  points(as.numeric(colnames(opt_vpa$N_mat))[years>=2010],colSums(opt_vpa$N_mat)[years>=2010]/1000,type="b",pch=21,col="red",bg="pink")
  axis(1)
legend("topright",legend = c("Simple Mod","Gulland Approach"),pch=c(21,21),pt.bg=c("lightblue","pink"),col=c("steelblue","red"),lty=1)

#Age3+ Abundance
plot(years,rowSums(N_array[,3:6,"0.42"])/1000, type = "b", pch = 21, col = "steelblue", bg = "lightblue",ylim=c(0,3000),ylab="Abund. 3+ (x1,000)",main="VPA Est. Age3+ YEP")
  points(as.numeric(colnames(opt_vpa$N_mat)),colSums(opt_vpa$N_mat[3:6,])/1000,type="b",pch=21,col="red",bg="pink")
legend("topright",legend = c("Simple Mod","Gulland Approach"),pch=c(21,21),pt.bg=c("lightblue","pink"),col=c("steelblue","red"),lty=1)

#Age3+ Abundance - 2010 - present
plot(years[years>=2010],rowSums(N_array[years>=2010,3:6,"0.42"])/1000, type = "b", pch = 21, col = "steelblue", bg = "lightblue",ylim=c(0,600),xlab="years",ylab="Abund. 3+ (x1,000)",main="VPA Est. Age3+ YEP")
  points(as.numeric(colnames(opt_vpa$N_mat))[years>=2010],colSums(opt_vpa$N_mat[3:6,])[years>=2010]/1000,type="b",pch=21,col="red",bg="pink")
legend("topright",legend = c("Simple Mod","Gulland Approach"),pch=c(21,21),pt.bg=c("lightblue","pink"),col=c("steelblue","red"),lty=1)

#Age 1 Abundance
plot(years,N_array[,1,"0.42"]/1000, type = "b", pch = 21, col = "steelblue", bg = "lightblue",ylim=c(0,2500),ylab="Abund. Age1 (x1,000)",main="VPA Est. Age1 YEP")
  points(as.numeric(colnames(opt_vpa$N_mat)),opt_vpa$N_mat[1,]/1000,type="b",pch=21,col="red",bg="pink")
legend("topright",legend = c("Simple Mod","Gulland Approach"),pch=c(21,21),pt.bg=c("lightblue","pink"),col=c("steelblue","red"),lty=1)

#Age3+ Abundance - 2010 - present
plot(years,N_array[,1,"0.42"]/1000, type = "b", pch = 21, col = "steelblue", bg = "lightblue",ylim=c(0,150),ylab="Abund. Age1 (x1,000)",main="VPA Est. Age1 YEP",xlim=c(2010,2025))
  points(as.numeric(colnames(opt_vpa$N_mat)),opt_vpa$N_mat[1,]/1000,type="b",pch=21,col="red",bg="pink")
legend("topright",legend = c("Simple Mod","Gulland Approach"),pch=c(21,21),pt.bg=c("lightblue","pink"),col=c("steelblue","red"),lty=1)

#Plot Age5/6 F
# plot(NA,NA, type = "n", pch = 21, col = "red", bg = "pink",xlim=c(min(years),max(years)),ylim=c(0,2),ylab="Z",main="Age 5 and 6",las=2,xaxt="n")
#   axis(1)
#   polygon(c(years,rev(years)),c(opt_vpa$F_mat[6,]+M,rep(0,length(years))),col = "pink",border=NA)
#   polygon(c(years,rev(years)),c(rep(M,length(years)),rep(0,length(years))),col = "darkgrey",border = NA)
# legend("topright",legend = c("F","M"),fill=c("pink","darkgrey"))

#Plot as a barplot instaed
x<-barplot(rbind(M,opt_vpa$F_mat[6,]),col=c("darkgrey","pink"),xaxt="n",las=2,ylim=c(0,2),ylab="Z")
  axis(1,at = x[years%in%seq(1990,2020,5)],labels = years[years%in%seq(1990,2020,5)])
legend("topright",legend = c("F","M"),fill=c("pink","darkgrey"))


#How does the age3+ plot compare to the nets CPUE?
#Age3+ Abundance
par(mar=c(4.5,6,2,6))
plot(years,rowSums(N_array[,3:6,"0.42"])/1000, type = "b", pch = 21, col = "steelblue", bg = "lightblue",ylim=c(0,3000),ylab="Abund. 3+ (x1,000)",main="VPA Est. Age3+ YEP",las=2,xaxt="n")
  points(as.numeric(colnames(opt_vpa$N_mat)),colSums(opt_vpa$N_mat[3:6,])/1000,type="b",pch=21,col="red",bg="pink")
  axis(1)
  
  par(new=TRUE)
  plot(years,rowSums(GN_mat[-nrow(GN_mat),3:6]), type = "b", pch = 22, col = "black", bg = "darkgrey",ylim=c(0,150),ylab="",main="",yaxt="n",xaxt="n")
  mtext("Nearshore GN CPUE (#Age3+/Net)",side = 4,line = 3.5)
  axis(4,las=2)
legend("topright",legend = c("Simple Mod","Gulland Approach","GN CPUE"),pch=c(21,21,22),pt.bg=c("lightblue","pink","darkgrey"),col=c("steelblue","red","black"),lty=1)


par(mar=c(4.5,6,2,6))
plot(years,rowSums(N_array[,3:6,"0.42"])/1000, type = "b", pch = 21, col = "steelblue", bg = "lightblue",ylim=c(0,3000),ylab="Abund. 3+ (x1,000)",main="VPA Est. Age3+ YEP",las=2,xaxt="n")
  points(as.numeric(colnames(opt_vpa$N_mat)),colSums(opt_vpa$N_mat[3:6,])/1000,type="b",pch=21,col="red",bg="pink")
  axis(1)
  
  par(new=TRUE)
  plot(years[years>=1999],rowSums(offGN_mat[-nrow(offGN_mat),3:6]), type = "b", pch = 22, col = "black", bg = "darkgrey",xlim=c(min(years),max(years)),ylim=c(0,100),ylab="",main="",yaxt="n",xaxt="n")
  mtext("Offshore GN CPUE (#Age3+/Net)",side = 4,line = 3.5)
  axis(4,las=2)
legend("topright",legend = c("Simple Mod","Gulland Approach","GN CPUE"),pch=c(21,21,22),pt.bg=c("lightblue","pink","darkgrey"),col=c("steelblue","red","black"),lty=1)

#Mismatch at peaks before and at 2000 - Switch to otoliths?


#######################################
#Sensitivity plots for the GUlland Mod#
#######################################
#How does Age3+ abundance move with different values of Fterm?

F_termz_abund<-matrix(NA,nrow=6,ncol = length(years))
rownames(F_termz_abund)<-c(0.05,0.25,opt_val$par,0.75,1,1.5)
colnames(F_termz_abund)<-years
F_termz_Fval<-F_termz_abund


for(z in 1:nrow(F_termz_abund)){
  vpa_valz<-Gulland_FUN(catch = t(ca),F_term = as.numeric(rownames(F_termz_abund))[z],opt = FALSE)  
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

retro_abund<-matrix(NA,nrow=6,ncol = length(years))
rownames(retro_abund)<-seq(0,5)
colnames(retro_abund)<-years
retro_F<-retro_abund

for(z in 1:nrow(retro_abund)){
  sim_years<-years[1:(length(years)-as.numeric(rownames(retro_abund)[z]))]
  opt_val<-nlminb(start = list(F_term = 2), #If initital value is <0.4, optimization will push terminal F to 0. 
        objective = Gulland_FUN, 
        catch = t(ca)[,as.character(sim_years)], 
        opt = TRUE,
        control = list(trace = 1, eval.max = 1000, iter.max = 1000),lower = 0.1,upper = Inf)
  opt_vpa<-Gulland_FUN(catch = t(ca)[,as.character(sim_years)],F_term = opt_val$par,opt = FALSE)
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

#Write the VPA abundance as a .rds
opt_val<-nlminb(start = list(F_term = 2), #If initital value is <0.4, optimization will push terminal F to 0. 
        objective = Gulland_FUN, 
        catch = t(ca), 
        opt = TRUE,
        control = list(trace = 1, eval.max = 1000, iter.max = 1000),lower = 0.1,upper = Inf)
opt_vpa<-Gulland_FUN(catch = t(ca),F_term = opt_val$par,opt = FALSE)
#saveRDS(object = opt_vpa,file = paste0(dat_path,"opt_vpa.rds"))


  

  
  




