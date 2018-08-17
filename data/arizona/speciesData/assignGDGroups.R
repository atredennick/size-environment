rm(list=ls(all=TRUE)); graphics.off();
require(mgcv); require(lme4);

state<-"Arizona"
i<-2
dataDir<-paste0("C:/repos/driversdata/Data/",state,"/speciesData")

sppList=sort(c("BOER","BORO"))

# ######Survival
# dat<-read.csv(paste0(dataDir,"/",sppList[i],"/survD_reduced.csv"),header=TRUE)
# 
# levels(dat$quad)
# dat$Group<-substr(dat$quad, 1, 1)
# 
# write.csv(dat,paste0(dataDir,"/",sppList[i],"/survD_reduced.csv"))
# 
# 
# #####Growth
# dat<-read.csv(paste0(dataDir,"/",sppList[i],"/growDnoNA_reduced.csv"),header=TRUE)
# 
# levels(dat$quad)
# dat$Group<-substr(dat$quad, 1, 1)
# 
# write.csv(dat,paste0(dataDir,"/",sppList[i],"/growDnoNA_reduced.csv"))



#####recruits
dat<-read.csv(paste0(dataDir,"/",sppList[i],"/recArea_reduced.csv"),header=TRUE)

levels(dat$quad)
dat$Group<-substr(dat$quad, 1, 1)

write.csv(dat,paste0(dataDir,"/",sppList[i],"/recArea_reduced.csv"))
