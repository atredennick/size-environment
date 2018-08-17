fetchDemoData<-function(doSpp,state,surv=TRUE,dataDir,sppList){
  
  # For exponential W when spline kernel could not be fitted (2 Kansas species) 
     alphaList <- c(0.3558,.01137);
     midRings <- c(seq(1,19,2),seq(22.5,47.5,5),seq(55,145,10));
  
  ##############################################################
  # Load demographic data
  ##############################################################
  
  #------------------------------ load data for growth 
  if(surv==FALSE){ 
  
  # load growth data
  nonCompLength.g=5; # Number of columns in SppData that are not measures of competitors
  
  if(state=="idaho"| state=="montana"){
       groDfile=paste(dataDir,"/speciesData/",doSpp,"/growDnoNA.csv",sep="")
    }else{
       groDfile=paste(dataDir,"/speciesData/",doSpp,"/growDnoNA_reduced.csv",sep="") 
  }
  
  groD=read.csv(file=groDfile)
  D=groD[groD$allEdge==0,]; 
  D$logarea.t0=log(D$area.t0)
  D$logarea.t1=log(D$area.t1)
  D$quad=as.character(D$quad)
  D=D[order(D$X),]
  sort(table(D$area.t0),decreasing=TRUE)[1:3]
  
  ##########################################################
  # Read in data on neighbors 
  ##########################################################
  ringD <- read.csv(paste(dataDir,"/speciesdata/",doSpp,"/",doSpp,"_nbhood_rings.csv",sep=""))
  ringD$year<-as.factor(ringD$year)
  
  # merge D with ringD (D contains fewer rows)
  D<-merge(D,ringD,by.x=c("quad","year","trackID"),by.y=c("quad","year","genetID"))
  D=D[order(D$X),]
  rm(ringD)
  row.names(D) <- NULL  
  
  #-- ANNULUS DATA: all species 
  allCols= which(substr(names(D),1,4)%in%sppList); 
  sppData<-data.frame(logarea.t1=D$logarea.t1,quad=as.factor(D$quad),year=D$year,ID=D$trackID,age=D$age,
            logarea.t0=D$logarea.t0,Group=as.factor(D$Group),as.matrix(D[,allCols]))
  
  #----------------------- or, load data for survival 
  }else{
    
    #--------- Idaho or Montana, all the data can be used. 
    nonCompLength.s=5 #Number of columns in SppData that are not measures of competitors 
    if(state=="idaho"| state=="montana"){
            survDfile=paste(dataDir,"/speciesData/",doSpp,"/survD.csv",sep="")
    #--------- Otherwise use a reduced data set, edited by Adler Postdoc Chen.  
      }else{   
            survDfile=paste(dataDir,"/speciesData/",doSpp,"/survD_reduced.csv",sep="") 
    }
    
    survD=read.csv(file=survDfile)
    D=survD[survD$allEdge==0,];
    D$year=D$year
    D$logarea=log(D$area)
    D$quad=as.character(D$quad)
    D1=D=D[order(D$X),];
    
  ##########################################################
  # Read in data on neighbors 
  ##########################################################
    ringD <- read.csv(paste(dataDir,"/speciesdata/",doSpp,"/",doSpp,"_nbhood_rings.csv",sep=""))
    ringD$year<-ringD$year
    
    # merge D with ringD (D contains fewer rows)
    D<-merge(D,ringD,by.x=c("quad","year","trackID"),by.y=c("quad","year","genetID"))
    D=D[order(D$X),]
    rm(ringD)
    row.names(D) <- NULL  
    
    # -- ANNULUS DATA: all species  
    allCols= which(substr(names(D),1,4)%in%sppList); 
    if(state=="idaho"){
      sppData<-data.frame(survives=D$survives,age=D$age,ID=D$trackID, year=D$year, logarea=D$logarea, 
                          Group=as.factor(D$Group), seedling=D$seedling, quad=D$quad, as.matrix(D[,allCols]));
    }else{
      sppData<-data.frame(survives=D$survives,age=D$age,ID=D$trackID, year=D$year, logarea=D$logarea, 
                          Group=as.factor(D$Group), quad=D$quad, as.matrix(D[,allCols]));
    }
    
    }
    # ------------------ end the load data for survival
  
#--------------------------- end the demographic import 
  
#---------------------------- get the distance weights for intraspecific competition   
   dists <- read.csv(paste(dataDir,"/speciesdata/",state,"DistanceWeights.csv",sep="")); 
   if(state=="kansas" & doSpp=="BOHI" | state=="kansas" & doSpp=="SCSC"){
    sppNum <- if(doSpp=="BOHI"){1}else{2} 
    alpha <- alphaList[sppNum]; 
    dist_wts <- exp(-alpha*midRings); 
    dists[,doSpp] <- dist_wts;
  }else{
    dist_wts<- dists[,doSpp]; 
  }

# -------------------------- calculate W's using "response model": same kernel for all neighbors 
  nSpp = length(sppList); 
  W.resp = matrix(NA,nrow(sppData),nSpp); 
  for(j in 1:nSpp) { 
    sppCols = which(substr(names(sppData),1,4)==sppList[j]); 
    C <- data.matrix(sppData[,sppCols]) #matrix of areas in the annuli 
    W.resp[,j] <- C%*%dist_wts; 
  }
  
# -------------------------- calculate W's using "effects model": each neighbor has its own kernel    
  W.eff = matrix(NA,nrow(sppData),nSpp); 
  for(j in 1:nSpp) { 
    sppCols = which(substr(names(sppData),1,4)==sppList[j]); 
    C <- data.matrix(sppData[,sppCols]) #matrix of areas in the annuli 
    W.eff[,j] <- C%*%dists[,j]; 
  } 
  
# ---------------------- Assemble the data frame for growth modeling   
  if(surv==FALSE){
  dataG<-data.frame(sppData[,c("logarea.t1","year","ID","age","logarea.t0","Group","quad")], W.resp, W.eff)
  dataG$year<-as.numeric(paste(19,dataG$year,sep=""))
  
  dataG$Group <- as.factor(dataG$Group);
  data<-dataG
  
# ----------------------- Or, the data frame for survival modeling  
  }else{
    if(state=="idaho"){
      dataS<-data.frame(sppData[,c("survives","year","ID","logarea","age","Group","quad","seedling")], W.resp, W.eff)
    }else{
      dataS<-data.frame(sppData[,c("survives","year","ID","logarea","age","Group","quad")], W.resp, W.eff)
    }
    
    dataS$year<-as.numeric(paste(19,dataS$year,sep=""))
    dataS$Group <- as.factor(dataS$Group);
    dataS$logarea.t0 <- dataS$logarea; 
    data<-dataS
  }
  
  Wcols = which(substr(names(data),1,1)=="X");
  nspp <- length(sppList)
  names(data)[Wcols] <- c( paste0("W.resp",as.character(c(1:nspp))), paste0("W.eff",as.character(c(1:nspp))) ) 
  
  data$uniqID = paste(data$quad,data$ID,sep="_"); 
  
  return(data)
}

