fetchDemoData<-function(doSpp,state,grow=FALSE,surv=FALSE,recr=FALSE,dataDir){
  
  # Function to load and format demographic data from five datasets.
  # Written by Brittany Teller, modified by Andrew Tredennick

  #Out-dated W approach
  if(state=="Kansas" & doSpp=="BOHI"){
    alphaList <- c(0.3558,.01137);
    midRings <- c(seq(1,19,2),seq(22.5,47.5,5),seq(55,145,10));}
  
  ##############################################################
  #load demographic data
  ##############################################################
  if(grow==TRUE){
    # load growth data
    nonCompLength.g=5; # Number of columns in SppData that are not measures of competitors
    if(state=="Idaho"| state=="Montana"){groDfile=paste(dataDir,"/speciesdata/",doSpp,"/growDnoNA.csv",sep="")}else{
      groDfile=paste(dataDir,"/speciesdata/",doSpp,"/growDnoNA_reduced.csv",sep="") 
    }
    groD=read.csv(file=groDfile)
    D=groD[groD$allEdge==0,]; 
    #D$year=as.factor(D$year)
    D$logarea.t0=log(D$area.t0)
    D$logarea.t1=log(D$area.t1)
    D$quad=as.character(D$quad)
    D=D[order(D$X),]
    
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
    
    ## pull out annulus data for the focal species  
    sppCols=which(substr(names(D),1,4)==doSpp); 
    sppData<-data.frame(logarea.t1=D$logarea.t1,quad=as.factor(D$quad),year=D$year,ID=D$trackID,age=D$age,logarea.t0=D$logarea.t0,Group=as.factor(D$Group),as.matrix(D[,sppCols]))
    
  } #import growth data
  
  if(surv==TRUE){
    
    nonCompLength.s=5 #Number of columns in SppData that are not measures of competitors 
    if(state=="Idaho"| state=="Montana"){survDfile=paste(dataDir,"/speciesdata/",doSpp,"/survD.csv",sep="")}else{
      survDfile=paste(dataDir,"/speciesdata/",doSpp,"/survD_reduced.csv",sep="") 
    }
    survD=read.csv(file=survDfile)
    D=survD[survD$allEdge==0,];
    D$year=D$year
    D$logarea=log(D$area)
    D$quad=as.character(D$quad)
    D<-D[order(D$X),];
    
    #########read data on neighbors
    ringD <- read.csv(paste(dataDir,"/speciesdata/",doSpp,"/",doSpp,"_nbhood_rings.csv",sep=""))
    ringD$year<-ringD$year
    
    # merge D with ringD (D contains fewer rows)
    D<-merge(D,ringD,by.x=c("quad","year","trackID"),by.y=c("quad","year","genetID"))
    D=D[order(D$X),]
    rm(ringD)
    row.names(D) <- NULL  
    
    ## pull out annulus data for the focal species  
    sppCols=which(substr(names(D),1,4)==doSpp); 
    sppData<-data.frame(survives=D$survives,age=D$age,ID=D$trackID, year=D$year, logarea=D$logarea, Group=as.factor(D$Group), quad=D$quad, as.matrix(D[,sppCols]));

  } #Import survival data
  
  if(recr==TRUE){
    
    # get recruitment data
    if(state=="Idaho"| state=="Montana"){infile1=paste(dataDir,"/speciesdata/",doSpp,"/recArea.csv",sep="")}else{
      infile1=paste(dataDir,"/speciesdata/",doSpp,"/recArea_reduced.csv",sep="") 
    }
    tmpD=read.csv(infile1)
    tmpD=tmpD[,c("quad","year","NRquad","totParea","Group","recArea")]
    names(tmpD)[3]="R"
    names(tmpD)[4]="cov"
    D=tmpD
    D[is.na(D)]=0  # replace missing values 
    
    # calculate mean cover by group and year
    tmpD1=D[,c("quad","year","Group","cov")]
    tmpD1=aggregate(tmpD1$cov,by=list("year"=tmpD1$year,
                                    "Group"=tmpD1$Group),FUN=mean)
    names(tmpD1)[3]="Gcov"
    
    D=merge(D,tmpD1,all.x=T)
    D$y=D$R; 
    D$parents1=D$cov; D$parents2=D$Gcov; 
    D$fyear = factor(D$year); 
    intraD<-D
    
  } #Import recruitment data
  
  if(recr==FALSE){ #if this is growth or survival, compute W otherwise, don't.
    intraD<-sppData  # focal spp and intraspecific neighbors
    intraD<-data.frame(intraD)  
    
    if(state=="Kansas" & doSpp=="BOHI" | state=="Kansas" & doSpp=="SCSC"){
      sppNum <- if(doSpp=="BOHI"){1}else{2} 
      alpha <- alphaList[sppNum]; 
      dist_wts <- exp(-alpha*midRings); 
      
    }else{
      dists <- read.csv(paste(dataDir,"/speciesdata/",state,"DistanceWeights.csv",sep="")); 
      dist_wts<- dists[,paste0(doSpp)]
    }
    
    C <- data.matrix(intraD[,grep(paste(doSpp),names(intraD),value=T)]) #matrix of conspecific areas in the annuli 
    W <- C%*%dist_wts; 
    
    if(grow==TRUE){
      datag<-data.frame(intraD[,c("logarea.t1","year","ID","age","logarea.t0","Group","quad")], W)
      datag$year<-as.numeric(paste(19,datag$year,sep=""))
      datag$Group <- as.factor(datag$Group);
      data<-datag
    }
    if(surv==TRUE){
      
      datas<-data.frame(intraD[,c("survives","year","ID","logarea","age","Group","quad")], W)
      datas$year<-as.numeric(paste(19,datas$year,sep=""))
      datas$Group <- as.factor(datas$Group);
      data<-datas
    }#end survival T/F
    
  }# W import for surv and grow T/F merge
  
  if(recr==TRUE){ #for recruitment, compute true parents
    datar<-data.frame(intraD)
    datar$year<-as.numeric(paste(19,datar$year,sep=""))
    datar$Group <- as.factor(datar$Group);
    
    data<-datar
  } # and merge
  
  return(list(data=data))
}


  
  