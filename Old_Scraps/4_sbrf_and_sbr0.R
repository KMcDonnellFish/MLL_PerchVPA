#calculate lx

# life 1 history vectors
vul <- c(0.05, 0.174, 0.571, 0.899, 1.000, 0.6) # vulnerability at age
wa <- c(0.05, 0.104, 0.145, 0.186, 0.212, 0.253) # weight at age
mat <- c(0.05, 0.20, 0.82, 0.99, 1.00, 1.00) # maturity at age
M <- 0.4 # instantaneous natural mortality
F <- 0.1 # instantaneous fishing mortality
agez<-length(vul)

for(i in 1:agez){
  if(i==1) {
    lx<-numeric(agez)
    lx[1]<-1
    } else {
      lx[i]<-lx[i-1]*exp(-(F*vul[i-1]+M))
      }
    if(i==agez){ #Correct for plus group.
      lx[i]<-lx[i-1]/(1-exp(-(F*vul[i]+M)))
      }
  }

#Calcualte sbrf
sum(lx*mat*wa)


# now do it for Mille Lacs YEP

#Set Paths
dat_path<-"C:/Users/EU01237640/OneDrive - State of Minnesota - MN365/Mille_Lacs/YEP/PerchVPA/"

#I only have length/weight data from 1999 to present.  THey used scales prior to 1999 and I do not have those data.

#Load data
vpa_dat<-readRDS(file = paste0(dat_path,"opt_vpa.rds"))
LW_dat<-readRDS(file = paste0(dat_path,"LW_byage.rds"))
mat_dat<-readRDS(file = paste0(dat_path,"mat_byage.rds"))

LE_vul<-rbind(Unit_1=c(0.212,0.675,1,0.996,0.295),
                Unit_2=c(0.09,0.432,0.812,1,1),
                Unit_3=c(0.03,0.243,0.611,0.842,1),
                Unit_4=c(0.1,0.445,0.885,1,0.662))
colnames(LE_vul)<-seq(2,6)

#Assume Age 1 vul is 50% of 
#Check to see how much smaller age 1 fish are then age 2 ~50% smaller.
mean(LW_dat[,1,"length"]/LW_dat[,2,"length"])
LE_vul<-cbind("1"=LE_vul[,"2"]*0.5,LE_vul)
#Take mean across units
vul<-apply(LE_vul,MARGIN = 2,FUN = mean)
#Now standardize by age 5
vul<-vul/vul["5"]
wa<-LW_dat[,as.character(seq(1,6,1)),"weight"]
mat<-mat_dat[,as.character(seq(1,6,1))]

M <- 0.2
F <- t(vpa_dat$F_mat)

sbr_mat<-data.frame(sbrf=rep(0,length(1999:2024)),
                    sbr0=rep(0,length(1999:2024)),row.names = 1999:2024)

#Now loop through each year and calculate sbrf and sbr0
for(t in 1999:2024){
  
  agez<-length(vul)
  
  #sbrf
  for(i in 1:agez){
    if(i==1) {
      lx<-numeric(agez)
      lx[1]<-1
      } else {
        lx[i]<-lx[i-1]*exp(-(F[as.character(t),i-1]*vul[i-1]+M))
        }
      if(i==agez){ #Correct for plus group.
        lx[i]<-lx[i-1]/(1-exp(-(F[as.character(t),i]*vul[i]+M)))
        }
    }
  #Calculate sbrf (kg)
  sbr_mat[as.character(t),"sbrf"]<-sum(lx*wa[as.character(t),]/1000*mat[as.character(t),])
  
  #sbr0
  for(i in 1:agez){
    if(i==1) {
      l0<-numeric(agez)
      l0[1]<-1
      } else {
        l0[i]<-l0[i-1]*exp(-(M))
        }
      if(i==agez){ #Correct for plus group.
        l0[i]<-l0[i-1]/(1-exp(-(M)))
        }
    }
  #Calculate sbr0 (kg)
  sbr_mat[as.character(t),"sbr0"]<-sum(l0*wa[as.character(t),]/1000*mat[as.character(t),])
  }


# plot it out...
windows()
plot(as.numeric(rownames(sbr_mat)),sbr_mat[,"sbr0"],las=2,xaxt="n",xlab="Year",ylab="sbr",type="l",ylim=c(0,2),
  lwd=1.5,col="steelblue")
  points(as.numeric(rownames(sbr_mat)),sbr_mat[,"sbr0"],pch=21,col="black",bg="steelblue")  
  axis(1,at=seq(2000,2024,2))
  axis(1,at=seq(1999,2024,1),labels = NA)
  lines(as.numeric(rownames(sbr_mat)),sbr_mat[,"sbrf"],lwd=1.5,col="darkred")
  points(as.numeric(rownames(sbr_mat)),sbr_mat[,"sbrf"],pch=21,col="black",bg="darkred")  
legend("topleft",legend=c("sbr0","sbrf"),pch=c(21,21),col=c("black","black"),pt.bg=c("steelblue","darkred"))

#saveRDS(object = sbr_mat,file = paste0(dat_path,"sbr_mat.rds"))






