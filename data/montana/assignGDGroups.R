rm(list=ls(all=TRUE)); graphics.off();
require(mgcv); require(lme4);

state<-"Montana"
i<-4
dataDir<-paste0("C:/repos/driversdata/Data/",state,"/speciesData")
sppList=sort(c("BOGR","HECO","PASM","POSE"))

# dat<-read.csv(paste0(dataDir,"/",sppList[i],"/survD.csv"),header=TRUE)
# 
# levels(dat$quad)
# dat$Group<-substr(dat$quad, 1, 1)
# 
# write.csv(dat,paste0(dataDir,"/",sppList[i],"/survD.csv"))
# 
# 
# dat<-read.csv(paste0(dataDir,"/",sppList[i],"/growDnoNA.csv"),header=TRUE)
# 
# levels(dat$quad)
# dat$Group<-substr(dat$quad, 1, 1)
# 
# write.csv(dat,paste0(dataDir,"/",sppList[i],"/growDnoNA.csv"))

dat<-read.csv(paste0(dataDir,"/",sppList[i],"/recArea.csv"),header=TRUE)

levels(dat$quad)
dat$Group<-substr(dat$quad, 1, 1)

write.csv(dat,paste0(dataDir,"/",sppList[i],"/recArea.csv"))
