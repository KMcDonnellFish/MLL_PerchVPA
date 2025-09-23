######################
#Set Up data----
######################
require("xlsx")
require("RTMB")

base_path<-getwd() #May have to manually change this if you don't open as folder 
dat_path<-paste0(base_path,"/Data/")
code_path<-paste0(base_path,"/Scripts/")
out_path<-paste0(base_path,"/Output/")

#Pull data out of the .dat file from original GLIFWC scaa model
datfile<-readLines(con = paste0(dat_path,"Perch2025.dat"))
GN_bio_dater<-read.xlsx(file = paste0(dat_path,"ML YEP fall GN.xlsx"),sheetIndex = "Sheet1",colIndex = 2:9) #Grab most recent bio data from this terrible xlsx file

#Scrub Data----
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

#Extracted harvest data, GN data, and offshore GN data. I not going to use the GN data formally, but will use it to compare
#vpa estimates to survey catch rates
catch_mat
GN_mat
offGN_mat


#Create Data Summaries----
#find colors for ages
colz<-c("#fde725","#5ec962","#21918c","#3b528b","#d093e0ff","#440154")

catch_plot<-apply(X = catch_mat,MARGIN = 1,FUN = cumsum)
jpeg(filename = paste0(out_path,"MLL_perch_harvest.jpeg"),width = 9,height = 7,units = "in",res = 250)
  par(mar=c(4.5,7,2,2))
  plot(0,0,type="n",xlim=c(1985,2025),ylim=c(0,1200000),cex.axis=1.2,las=2,xaxt="n",xlab="Year",ylab="",cex.lab=1.2)
  axis(1,cex.axis=1.2)
  plot_yearz<-as.numeric(dimnames(catch_plot)[[2]])
  for(i in nrow(catch_plot):1){
    polygon(x = c(plot_yearz,rev(plot_yearz)),y = c(catch_plot[i,],rep(0,length(catch_plot[i,]))),col=colz[i],border=NA)
    }
  mtext(text = "# YEP Harvested",side = 2,line = 5.5,cex=1.2)
  legend("topright",legend=paste0("age",seq(1,6)),fill=colz,cex=1.2)
dev.off()

GN_plot<-apply(X = GN_mat,MARGIN = 1,FUN = cumsum)
jpeg(filename = paste0(out_path,"MLL_perch_GNSurvey.jpeg"),width = 9,height = 7,units = "in",res = 250)
  par(mar=c(4.5,7,2,2))
  plot(0,0,type="n",xlim=c(1985,2025),ylim=c(0,250),cex.axis=1.2,las=2,xaxt="n",xlab="Year",ylab="",cex.lab=1.2)
  axis(1,cex.axis=1.2)
  plot_yearz<-as.numeric(dimnames(GN_plot)[[2]])
  for(i in nrow(GN_plot):1){
    polygon(x = c(plot_yearz,rev(plot_yearz)),y = c(GN_plot[i,],rep(0,length(GN_plot[i,]))),col=colz[i],border=NA)
    }
  mtext(text = "CPUE Nearshore Fall Survey",side = 2,line = 5.5,cex=1.2)
  legend("topright",legend=paste0("age",seq(1,6)),fill=colz,cex=1.2)
dev.off()

offGN_plot<-apply(X = offGN_mat,MARGIN = 1,FUN = cumsum)
jpeg(filename = paste0(out_path,"MLL_perch_offshoreGNSurvey.jpeg"),width = 9,height = 7,units = "in",res = 250)
  par(mar=c(4.5,7,2,2))
  plot(0,0,type="n",xlim=c(1985,2025),ylim=c(0,250),cex.axis=1.2,las=2,xaxt="n",xlab="Year",ylab="",cex.lab=1.2)
  axis(1,cex.axis=1.2)
  plot_yearz<-as.numeric(dimnames(offGN_plot)[[2]])
  for(i in nrow(offGN_plot):1){
    polygon(x = c(plot_yearz,rev(plot_yearz)),y = c(offGN_plot[i,],rep(0,length(offGN_plot[i,]))),col=colz[i],border=NA)
    }
  mtext(text = "CPUE Offshore Fall Survey",side = 2,line = 5.5,cex=1.2)
  legend("topright",legend=paste0("age",seq(1,6)),fill=colz,cex=1.2)
dev.off()

#Write data as .rds files for next step
saveRDS(object = catch_mat,file = paste0(dat_path,"catch_mat.rds"))
saveRDS(object = GN_mat,file = paste0(dat_path,"GN_mat.rds"))
saveRDS(object = offGN_mat,file = paste0(dat_path,"offGN_mat.rds"))

 #Length-at-age modeling ----
vonB_dater<-GN_bio_dater[,c("year","length..mm.","weight..gms.","age")]
colnames(vonB_dater)<-c("year","length","weight","age")
vonB_dater<-vonB_dater[!is.na(vonB_dater[,"age"]),]
vonB_dater$cohort<-vonB_dater[,"year"]-vonB_dater[,"age"]
#Remove any NA values
vonB_dater<-vonB_dater[!is.na(vonB_dater[,"length"]) & !is.na(vonB_dater[,"weight"]) & !is.na(vonB_dater[,"age"]) & !is.na(vonB_dater[,"cohort"]),]
#Make year and cohort factors
vonB_dater$year<-as.factor(vonB_dater[,"year"])
vonB_dater$cohort<-as.factor(vonB_dater[,"cohort"])

#Set up holding array
LW_array<-array(NA,dim = c(length(unique(vonB_dater[,"year"])),length(unique(vonB_dater[,"age"])),2),
              dimnames = list(year=unique(vonB_dater[,"year"]),
                              age=sort(unique(vonB_dater[,"age"])),
                              variable=c("length","weight")))

#Load Von B FUN
source(paste0(code_path,"1b_VonB.R"))

for(i in dimnames(LW_array)[[1]]){
  vonB_dater_singleyear<-vonB_dater[as.character(vonB_dater[,"year"])==i,]
  
  if(nrow(vonB_dater_singleyear)>10){
    params<-list(Linf=350,
                  k=0.2,
                  t0=-2,
                  L_pred_sd=10
                  )
    VonBFUN(params = params)
    obj<-MakeADFun(func = VonBFUN,parameters = params)
    
    bounds_mat<-data.frame(param_name=names(x = params),lower=-Inf,upper=Inf)
      bounds_mat[bounds_mat[,"param_name"]=="t0",c("lower","upper")]<-c(-Inf,-0.1)
    opt_lower<-bounds_mat[match(x = names(obj$par), table=bounds_mat[,"param_name"]),"lower"]
    opt_upper<-bounds_mat[match(x = names(obj$par), table=bounds_mat[,"param_name"]),"upper"]  
  
    opt<-nlminb(start = obj$par,objective = obj$fn,gradient = obj$gr,
                  control = list(eval.max=10000,iter.max=10000),lower = opt_lower,upper = opt_upper)
    LW_array[i,,"length"] <- opt$par["Linf"]*(1-exp(-opt$par["k"]*(as.numeric(dimnames(LW_array)[[2]])+0.4164 - opt$par["t0"]))) #Note this offset is to account for the ~5 months from beginning of year to time of capture during fall survey (~5 months from hatch).
    } 
  }

LW_array[,,"length"]


#Weight at age---- 
source(paste0(code_path,"1c_LW_reg.R"))
for(i in dimnames(LW_array)[[1]]){
  vonB_dater_singleyear<-vonB_dater[as.character(vonB_dater[,"year"])==i,]
  
  if(nrow(vonB_dater_singleyear)>10){
    params<-list(log_a=-10,
                  b=0,
                  w_sd=1)
    LW_reg(params = params)
    obj<-MakeADFun(func = LW_reg,parameters = params)
    opt<-nlminb(start = obj$par,objective = obj$fn,gradient = obj$gr,
                  control = list(eval.max=10000,iter.max=10000))
    LW_array[i,,"weight"] <- exp(opt$par["log_a"] + opt$par["b"]*log(LW_array[i,,"length"]))
    } 
  }

#Save all 
saveRDS(object = LW_array,file = paste0(dat_path,"LW_byage.rds"))

#Maturity at age ----
GN_sex_dater<-GN_bio_dater[!is.na(GN_bio_dater[,"sex"]),]
GN_sex_dater$mature<-0
unique(GN_sex_dater[,"sex"]) #Remove unknowns - assuming a "male uknown" is not mature.
GN_sex_dater<-GN_sex_dater[GN_sex_dater[,"sex"]%in%c("FM","MM","FI","MI","FU","MU","UI"),]

GN_sex_dater[GN_sex_dater[,"sex"]=="FM" | GN_sex_dater[,"sex"]=="MM","mature"]<-1

maturity <- array(NA,dim = dim(LW_array)[1:2],dimnames = dimnames(LW_array)[1:2])

for(i in dimnames(LW_array)[[1]]){
  glm_fit<-glm(mature~length..mm.,family=binomial,data= GN_sex_dater[GN_sex_dater[,"year"]==i,])
  maturity[i,] <- round(x = predict(object = glm_fit,newdata = data.frame(length..mm.= LW_array[i,,"length"]),type="response"),digits = 3)
  }

saveRDS(object = maturity,file = paste0(dat_path,"mat_byage.rds"))

#Lorenzen M calculation----
#Pool length data to get VonB parameters to estiamte Lorenzen M
vonB_dater[vonB_dater[,"age"]>=6,"age"]<-6
vonB_dater_singleyear<-vonB_dater
  
params<-list(Linf=350,
            k=0.2,
            t0=-2,
            L_pred_sd=10)
VonBFUN(params = params)
obj<-MakeADFun(func = VonBFUN,parameters = params)

bounds_mat<-data.frame(param_name=names(x = params),lower=-Inf,upper=Inf)
  bounds_mat[bounds_mat[,"param_name"]=="t0",c("lower","upper")]<-c(-Inf,-0.1)
opt_lower<-bounds_mat[match(x = names(obj$par), table=bounds_mat[,"param_name"]),"lower"]
opt_upper<-bounds_mat[match(x = names(obj$par), table=bounds_mat[,"param_name"]),"upper"]  

opt<-nlminb(start = obj$par,objective = obj$fn,gradient = obj$gr,
              control = list(eval.max=10000,iter.max=10000),lower = opt_lower,upper = opt_upper)

opt

pred_lengths<-opt$par["Linf"]*(1-exp(-opt$par["k"]*((1:6)+0.4164 - opt$par["t0"])))
#Use equation from Lorenzen 2022, Table 3: "Best Model"
#lnM = a + b*ln(L/Linf) + c*ln(k)
a <- 0.28
b <- -1.30
c <- 1.08

LorenzenM<-exp(a + b*log(pred_lengths/opt$par["Linf"]) + c*log(opt$par["k"]))
saveRDS(object = LorenzenM,file = paste0(dat_path,"LorenzenM.rds"))




