#########################################################################
### Exploratory data analysis to compare seedlings and others
### Seedlings are "likely seedlings" as defined by PBA 
### Code written by SPE and modified by ATT
#########################################################################
rm(list=ls(all.names=TRUE)); 
graphics.off();
require(mgcv); 
require(lme4); 
require(gamm4); 

## Provide directory information
root<-ifelse(.Platform$OS.type=="windows","c:/repos","~/repos"); # modify as needed
setwd(paste0(root,"/drivers/empirical/size_by_year_IPM/seedling_analysis/")); 

source("fetchDemoData-allWs.R") 
source("trimQuadrats.R"); 

state_species <- data.frame(state_name = c(rep("idaho",4),
                                      rep("montana",4),
                                      rep("arizona",2),
                                      rep("newmexico",2),
                                      rep("kansas",3)),
                            species = c(sort(c("ARTR","HECO","POSE","PSSP")),
                                        sort(c("BOGR","HECO","PASM","POSE")),
                                        sort(c("BOER","BORO")),
                                        sort(c("BOER","SPFL")),
                                        sort(c("ANGE","BOCU","BOHI"))))
state_species$species <- as.character(state_species$species)
all_states <- (unique(state_species$state))

for(state in all_states){
  sppList <- subset(state_species, state_name==state)$species
  ########################### Survival ##############################################
  Sfits=list(length(sppList));
  for(species in 1:length(sppList)){
    ## Import Demographic Data
    doSpp=sppList[species]; 
    dataDir<-paste0(root,"/driversdata/Data/",state, sep="")
    data<-fetchDemoData(doSpp=sppList[species], state=state, surv=TRUE, dataDir=dataDir,sppList=sppList)
    seedling_size <- round(as.numeric(names(sort(table(exp(data$logarea)),decreasing=TRUE)[1])),5)
    if(seedling_size > 0.26) seedling_size <- 0.26
    out=trimQuadrats(data, seedling_size); 
    data=out$data; 
    data$year<-as.factor(data$year)
    
    likely.seedling = which((data$age==1)&(data$logarea.t0<=log(seedling_size))) #catches all seedlings 
    cat(doSpp,out$doubtful,out$seedling,length(likely.seedling),"\n")
    
    if(state!="idaho"){
      data$seedling <- 0
    }
    data$seedling[likely.seedling] <- 1; 
    data$seedling <- factor(data$seedling); 
    
    wcols <- which(substr(names(data),1,3)=="W.r"); 
    W = data[,wcols]; data$W=W[,species]; 
    
    if(length(which(data$seedling==1))>0){
      gamS <- gam(survives ~ logarea + W + Group + seedling + seedling:Group + s(year,bs="re") + s(year,by=logarea,bs="re") + s(year,by=W,bs="re"), 
                  data=data, family=binomial) 
    }else{
      gamS <- gam(survives ~ logarea + W + Group + s(year,bs="re") + s(year,by=logarea,bs="re") + s(year,by=W,bs="re"), 
                  data=data, family=binomial) 
    }
    
    Sfits[[species]]=summary(gamS)$p.table;
    cat(species,"\n"); 
  } # end survival species loop
  
  ########################### Growth ##############################################
  source("fetchDemoData-allWs.R"); source("trimQuadrats.R"); 
  
  Gfits=list(length(sppList));
  source("fetchDemoData-allWs.R") 
  for(species in 1:length(sppList)){
    
    ## Import Demographic Data
    doSpp=sppList[species]; 
    dataDir<-paste0(root,"/driversdata/Data/",state, sep="")
    data<-fetchDemoData(doSpp=sppList[species], state=state, surv=FALSE, dataDir=dataDir,sppList=sppList)
    seedling_size <- round(as.numeric(names(sort(table(exp(data$logarea.t0)),decreasing=TRUE)[1])),5)
    if(seedling_size > 0.26) seedling_size <- 0.26
    out=trimQuadrats(data, seedling_size); data=out$data; 
    data$year<-as.factor(data$year)
    
    likely.seedling = which((data$age==1)&(data$logarea.t0<=log(0.26))) #catches all seedlings 
    cat(doSpp,out$doubtful,out$seedling,length(likely.seedling),"\n")
    
    data$year<-as.factor(data$year)
    
    data$seedling <- rep(0,nrow(data)); data$seedling[likely.seedling] <- 1; 
    data$seedling <- factor(data$seedling); 
    
    wcols <- which(substr(names(data),1,3)=="W.r"); 
    W = data[,wcols]; data$W=W[,species]; 
    
    if(length(which(data$seedling==1))>0){
      gamG <- gam(logarea.t1 ~ logarea.t0 + W + Group + seedling + seedling:Group + s(year,bs="re") + s(year,by=logarea.t0,bs="re") + s(year,by=W,bs="re"), 
                  data=data) 
    }else{
      gamG <- gam(logarea.t1 ~ logarea.t0 + W + Group +  s(year,bs="re") + s(year,by=logarea.t0,bs="re") + s(year,by=W,bs="re"), 
                  data=data) 
    }

    
    Gfits[[species]]=summary(gamG)$p.table; 
    if(length(likely.seedling>1)){
      #hist(data$logarea.t1[likely.seedling],main=doSpp); 
    }
    cat(species,"\n"); 
  } # end growth species loop
  
  names(Sfits) <- names(Gfits) <- sppList
  out_list <- list(Sfits=Sfits,Gfits=Gfits)
  saveRDS(out_list, paste0("seedling_test_fits_",state,".RDS"))
} # end state loop





