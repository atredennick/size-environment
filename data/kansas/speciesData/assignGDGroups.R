rm(list=ls(all=TRUE)); graphics.off();
require(mgcv); require(lme4);

state<-"Kansas"
i<-4
dataDir<-paste0("C:/repos/driversdata/Data/",state,"/speciesData")

setwd("C:/repos/driversdata/data/kansas/speciesData")
info<-read.csv("quadrat_info.csv",header=T)
subinfo<-data.frame(quad=info$shapefiles, Group=info$group, Grazing=info$Grazing)

sppList=sort(c("ANGE","BOCU","BOHI","SCSC"))

####Survival
dat<-read.csv(paste0(dataDir,"/",sppList[i],"/survD_reduced.csv"),header=TRUE)

new<-merge(dat, subinfo, by="quad")
new<-new[,-2]

write.csv(new,paste0(dataDir,"/",sppList[i],"/survD_reduced.csv"))

#####Growth
dat<-read.csv(paste0(dataDir,"/",sppList[i],"/growDnoNA_reduced.csv"),header=TRUE)

new<-merge(dat, subinfo, by="quad")
new<-new[,-2]

write.csv(new,paste0(dataDir,"/",sppList[i],"/growDnoNA_reduced.csv"))


#####Recruitment
dat<-read.csv(paste0(dataDir,"/",sppList[i],"/recArea_reduced.csv"),header=TRUE)

new<-merge(dat, subinfo, by="quad")
new<-new[,-2]

write.csv(new,paste0(dataDir,"/",sppList[i],"/recArea_reduced.csv"))
