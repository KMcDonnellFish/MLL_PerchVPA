#This script runs the most recent iteration of GLIFWC's SCAA model for MLL YEP.
#These data and model were provided by Adam Ray so we can compare SCAA output with VPA modeling effort.

#Set Paths----
base_path<-getwd() #May have to manually change this if you don't open as folder 
admb_path<-paste0(base_path,"/Scripts/SCAA_admb/")
out_path<-out_path<-paste0(base_path,"/Output/")
dat_path<-paste0(base_path,"/Data/")

#Load packages----
require("R2admb")

#Compile and run model----
#Have to set working directory to the where the .dat and .tpl files are located (named the same)...
setwd(admb_path)
compile_admb(fn = "Perch2025",verbose = TRUE)
run_admb(fn = "Perch2025",verbose = TRUE)

#Looks like there is some bounding issues with a few parameters, should probably investigate this if we intend to make 
# any further inference on this model

#Grab model results----
read_rep(fn = "Perch2025") #This report file is unformatted. going to pull out stuff by lines instead

##Abundance----
rep_text<-readLines(con = "Perch2025.rep")
abun_index<-c(which(grepl(x = rep_text,pattern = "Predicted Abundance"))+1,
                which(grepl(x = rep_text,pattern = "Three Plus Abundance By Year"))-1)

est_abund<-rep_text[abun_index[1]:abun_index[2]]
scaa_abund<-simplify2array(lapply(X = strsplit(x = est_abund,split = " "),FUN = as.numeric))
#Remove NAs
scaa_abund<-scaa_abund[-1,]
dimnames(scaa_abund)<-list(seq(1,dim(scaa_abund)[1]),seq(1986,2024+1)) #MAGIC numbers here indicating years last year is a projection....

##F values ----
F_index<-c(which(grepl(x = rep_text,pattern = "Instantaneous F Matrix"))+1,
                which(grepl(x = rep_text,pattern = "Observed Total Kill By Year"))-1)
est_F<-rep_text[F_index[1]:F_index[2]]
scaa_F<-simplify2array(lapply(X = strsplit(x = est_F,split = " "),FUN = as.numeric))
scaa_F<-scaa_F[-1,] #Remove NAs
dimnames(scaa_F)<-list(seq(1,dim(scaa_F)[1]),seq(1986,2024)) #MAGIC numbers here indicating years

years<-as.numeric(dimnames(scaa_abund)[[2]])

#Plot Abundance and F----
##Age3 abundance----
jpeg(filename = paste0(out_path,"SCAA_age3_abund.jpeg"),width = 9,height = 7,units = "in",res = 250)
  par(mar=c(4.5,4.5,2,2))
  plot(years,colSums(scaa_abund[3:6,])/1000, type = "b", pch = 21, col = "steelblue", bg = "darkblue",ylim=c(0,2000),ylab="Abund. 3+ (x1,000)",las=2,xaxt="n")
    axis(1)
    # points(as.numeric(colnames(opt_vpa$N_mat)),colSums(opt_vpa$N_mat[3:6,])/1000,type="b",pch=21,col="red",bg="pink")
#   legend("topright",legend = c("Simple Mod","Gulland Approach"),pch=c(21,21),pt.bg=c("lightblue","pink"),col=c("steelblue","red"),lty=1)
dev.off()

##Compare scaa and vpa abundance----
#Bring in VPA data
opt_vpa<-readRDS(file = paste0(dat_path,"opt_vpa.rds"))

jpeg(filename = paste0(out_path,"SCAA_vs_vpa_age3_abund.jpeg"),width = 9,height = 7,units = "in",res = 250)
  par(mar=c(4.5,4.5,2,2))
  plot(years,colSums(scaa_abund[3:6,])/1000, type = "b", pch = 21, col = "steelblue", bg = "darkblue",ylim=c(0,8000),ylab="Abund. 3+ (x1,000)",las=2,xaxt="n")
    axis(1)
    points(as.numeric(colnames(opt_vpa$N_mat)),colSums(opt_vpa$N_mat[3:6,])/1000,type="b",pch=21,col="red",bg="pink")
  legend("topright",legend = c("SCAA","VPA"),pch=c(21,21),pt.bg=c("darkblue","pink"),col=c("steelblue","red"),lty=1)
dev.off()

##F by age----
#This model assumes M is 0.3 for all ages
jpeg(filename = paste0(out_path,"scaa_mort_by_age.jpeg"),width = 9,height = 7,units = "in",res = 250)
  par(mfcol=c(3,2))
  for(i in 1:6){
    x<-barplot(rbind(0.3,scaa_F[i,]),col=c("darkgrey","pink"),xaxt="n",las=2,ylim=c(0,2.5),ylab="Z")
      axis(1,at = x[years%in%seq(1990,2020,5)],labels = years[years%in%seq(1990,2020,5)])
      legend("topleft",legend = c("F","M"),fill=c("pink","darkgrey"),title = paste0("Age ",i))
    }
dev.off()

#Without messing too much with the SCAA model, I wonder how if the SCAA model would rescale with an updated Lorenzen M?
#Updating the M in the SCAA----
#M values are inthe data file

##Alter M values -----
#Import lorenzen M
Lorenzen_M<-readRDS(file = paste0(dat_path,"LorenzenM.rds"))

#Grab dat file - save a copy to restore at end.
dat_file_new<-dat_file_original<-readLines(con = "Perch2025.dat")

M_vals<-as.numeric(strsplit(x = dat_file_new[which(grepl(x = dat_file_new,"Natural mortalities"))+1],split = " ")[[1]])
#Check vals
M_vals
Lorenzen_M

dat_file_new[which(grepl(x = dat_file_new,"Natural mortalities"))+1]<-paste(Lorenzen_M,collapse = " ")

#Overwrite datfile
writeLines(text = dat_file_new,con = "Perch2025.dat")

##compile/rerun model----
compile_admb(fn = "Perch2025",verbose = TRUE)
run_admb(fn = "Perch2025",verbose = TRUE)
#Got a lot more bounding errors on F, probably shouldn't read into this too much, but I really just want to see if it rescales

###Grab output----
rep_text<-readLines(con = "Perch2025.rep")
abun_index<-c(which(grepl(x = rep_text,pattern = "Predicted Abundance"))+1,
                which(grepl(x = rep_text,pattern = "Three Plus Abundance By Year"))-1)

est_abund<-rep_text[abun_index[1]:abun_index[2]]
scaa_abund<-simplify2array(lapply(X = strsplit(x = est_abund,split = " "),FUN = as.numeric))
#Remove NAs
scaa_abund<-scaa_abund[-1,]
dimnames(scaa_abund)<-list(seq(1,dim(scaa_abund)[1]),seq(1986,2024+1)) #MAGIC numbers here indicating years last year is a projection....

F_index<-c(which(grepl(x = rep_text,pattern = "Instantaneous F Matrix"))+1,
                which(grepl(x = rep_text,pattern = "Observed Total Kill By Year"))-1)
est_F<-rep_text[F_index[1]:F_index[2]]
scaa_F<-simplify2array(lapply(X = strsplit(x = est_F,split = " "),FUN = as.numeric))
scaa_F<-scaa_F[-1,] #Remove NAs
dimnames(scaa_F)<-list(seq(1,dim(scaa_F)[1]),seq(1986,2024)) #MAGIC numbers here indicating years

years<-as.numeric(dimnames(scaa_abund)[[2]])

###Make Plots----

jpeg(filename = paste0(out_path,"SCAA_vs_vpa_age3_abund_SCAAWithLorenzM.jpeg"),width = 9,height = 7,units = "in",res = 250)
  par(mar=c(4.5,4.5,2,2))
  plot(years,colSums(scaa_abund[3:6,])/1000, type = "b", pch = 21, col = "steelblue", bg = "darkblue",ylim=c(0,8000),ylab="Abund. 3+ (x1,000)",las=2,xaxt="n")
    axis(1)
    points(as.numeric(colnames(opt_vpa$N_mat)),colSums(opt_vpa$N_mat[3:6,])/1000,type="b",pch=21,col="red",bg="pink")
  legend("topright",legend = c("SCAA","VPA"),pch=c(21,21),pt.bg=c("darkblue","pink"),col=c("steelblue","red"),lty=1)
dev.off()

##F by age----
#This model assumes M is 0.3 for all ages
jpeg(filename = paste0(out_path,"scaa_mort_by_age_SCAAWithLorenzM.jpeg"),width = 9,height = 7,units = "in",res = 250)
  par(mfcol=c(3,2))
  for(i in 1:6){
    x<-barplot(rbind(Lorenzen_M[i],scaa_F[i,]),col=c("darkgrey","pink"),xaxt="n",las=2,ylim=c(0,2.5),ylab="Z")
      axis(1,at = x[years%in%seq(1990,2020,5)],labels = years[years%in%seq(1990,2020,5)])
      legend("topleft",legend = c("F","M"),fill=c("pink","darkgrey"),title = paste0("Age ",i))
    }
dev.off()

##Put the datfile back to the way it was and rerun----
writeLines(text = dat_file_original,con = "Perch2025.dat")
compile_admb(fn = "Perch2025",verbose = TRUE)
run_admb(fn = "Perch2025",verbose = TRUE)




