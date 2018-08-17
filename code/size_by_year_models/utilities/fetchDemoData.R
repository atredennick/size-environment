fetchDemoData<-function(doSpp,nUnits,WEEK,state,grow=FALSE,surv=FALSE,recr=FALSE,dataDir){

  #Out-dated W approach
  if(state=="Kansas" & doSpp=="BOHI"){
    alphaList <- c(0.3558,.01137);
    midRings <- c(seq(1,19,2),seq(22.5,47.5,5),seq(55,145,10));}
  
  #Simulation Parameters 
  monthStart=6;   
  
  #Climate Data
  DETREND=TRUE; 
  
  #Functions
  aggLags<-function(datC, fun, meas, seg){
    pars <- as.list(match.call()[-1])
    measure<-datC[,as.character(pars$meas)]
    segment<-datC[,as.character(pars$seg)]
    agg<-with(datC,ftable(tapply(measure, list(datC$year, segment), fun)))
    agg<-as.data.frame(agg)
    names(agg)<-c("year",pars$seg,paste(pars$meas, pars$fun,"Now", sep=""))
    
    return(list(agg,paste(pars$fun)))
  }
  
  appLags<-function(agDatC,nUnits,name, fun){
    storeNames<-names(agDatC)
    for (i in 1:nUnits){
      agDatC<-as.data.frame(agDatC)
      waved<-agDatC[[3]]
      newWave<-c(rep(NA, times=i),waved) 
      newThing<-newWave[1:nrow(agDatC)]
      agDatC<-cbind(agDatC, newThing)
    }
    names(agDatC)<-c(storeNames,paste(name,fun,1:nUnits, sep=""))
    return(agDatC)
  }
  
  #Functions revised for when climate depends on quad
  aggLagsQ<-function(datC, fun, meas, seg){
    pars <- as.list(match.call()[-1])
    measure<-datC[,as.character(pars$meas)]
    segment<-datC[,as.character(pars$seg)]
    agg<-with(datC,ftable(tapply(measure, list(datC$yearQ, segment), fun)))
    agg<-as.data.frame(agg)
    names(agg)<-c("yearQ",pars$seg,paste(pars$meas, pars$fun,"Now", sep=""))
    
    return(list(agg,paste(pars$fun)))
  }
  
  
  #Different states need to be formatted differently
  if(state=="Idaho" | state=="Montana"){  # daily climate data detrended & aggregated into lags
    
    #Load Climate Data
    climDfile=paste(dataDir,"/climateData/interpClim.csv",sep="")
    iClim=read.csv(file=climDfile)
    iClim<-iClim[,-1]
    
    #assign week of year
    iClim$week<-NA
    yearvec<-unique(iClim$year)
    
    for (i in 1:length(yearvec)){
      yeari<-yearvec[i]
      subs<-subset(iClim,iClim$year==yeari)
      week<-rep(1:53, each=7)
      iClim$week[iClim$year==yeari]<-week[1:nrow(subs)]  
    }
    
    ## Change units: Ppt in mm, temperatures in Celsius. 
    if(state=="Idaho"){
      iClim$Ppt<-iClim$Ppt/0.039370
      iClim$Tmax<-(5/9)*(iClim$Tmax-32);
      iClim$Tmin<-(5/9)*(iClim$Tmin-32);
      iClim$Tmid<-(iClim$Tmin+iClim$Tmax)/2; 
    }else{
      iClim$Tmid<-(iClim$Tmin+iClim$Tmax)/2
    }
    
    if(DETREND) {
      require(mgcv); #par(mfrow=c(2,2))
      iClim$DOY<- 30*(iClim$month-1)+iClim$day
      
      fit=gam(Ppt~s(DOY,bs="cc"),gamma=1.4,data=iClim,method="GCV.Cp")
      iClim$Ppt <- fit$residuals; #plot(fit,main="Precip");
      
      fit=gam(Tmax~s(DOY,bs="cc"),gamma=1.4,data=iClim,method="GCV.Cp")
      iClim$Tmax <- fit$residuals; #plot(fit,main="Tmax");
      
      fit=gam(Tmin~s(DOY,bs="cc"),gamma=1.4,data=iClim,method="GCV.Cp")
      iClim$Tmin <- fit$residuals; #plot(fit,main="Tmin");
      
      fit=gam(Tmid~s(DOY,bs="cc"),gamma=1.4,data=iClim,method="GCV.Cp")
      iClim$Tmid <- fit$residuals; #plot(fit,main="Tmax");
    }
    
    
    #######################################################################
    #Generate lags of climate signals
    #######################################################################
    
    #make monthly averages 
    tempAggWT<-aggLags(datC=iClim, fun=mean, meas=Tmid, seg=month)
    subAgg<-tempAggWT[[1]]
    subAgg$month<-as.numeric(as.character(subAgg$month))
    subAgg<-subAgg[with(subAgg,order(year,month)),]
    
    tempWaveWT<-appLags(agDatC=subAgg, nUnits=nUnits, name=paste("month","Tmid", sep=""), fun=tempAggWT[[2]])
    tempWaveWT<-na.omit(tempWaveWT)
    
    monthStart=monthStart
    subsWT<-subset(tempWaveWT,tempWaveWT$month==monthStart) #choose start month
    subsWT<-subsWT[,-which(names(subsWT)=="month")] #Get rid of "month" column
    subsWT$year<-as.numeric(as.character(subsWT$year)) #make "year" readable to merge
    
    #reshape to show temp over the last nUnits
    labs<-c("year","t.00", paste("t.0",1:9, sep=""),
            paste("t.",10:nUnits, sep=""))
    names(subsWT) <- labs
    
    tempAggWP<-aggLags(datC=iClim, fun=mean, meas=Ppt, seg=month)
    subAgg<-tempAggWP[[1]]
    subAgg$month<-as.numeric(as.character(subAgg$month))
    subAgg<-subAgg[with(subAgg,order(year,month)),]
    
    tempWaveWP<-appLags(agDatC=subAgg, nUnits=nUnits, name=paste("month","Ppt", sep=""), fun=tempAggWP[[2]])
    tempWaveWP<-na.omit(tempWaveWP)
    subsWP<-subset(tempWaveWP,tempWaveWP$month==monthStart) #choose a start month
    subsWP<-subsWP[,-which(names(subsWP)=="month")] #remove month column
    subsWP$year<-as.numeric(as.character(subsWP$year)) #make "year" readable to merge
    
    labs<-c("year","p.00",paste("p.0",1:9, sep=""),
            paste("p.",10:nUnits, sep=""))
    names(subsWP) <- labs
    
    climDat<-cbind(subsWP,subsWT)
    
  }
  
  if(state=="Kansas" | state=="NewMexico"){ #Monthly averages already available
    #Load Climate Data
    climDfile=paste(dataDir,"/climateData/monthClim.csv",sep="")
    iClim=read.csv(file=climDfile)
    iClim<-iClim[,-1]
    
    #Detrend by month
    if(DETREND) {
      require(mgcv); #par(mfrow=c(2,2))
      
      fit=gam(Ppt~s(month,bs="cc"),gamma=1.4,data=iClim,method="GCV.Cp")
      iClim$Ppt <- fit$residuals; #plot(fit,main="Precip");
      
      fit=gam(Tmid~s(month,bs="cc"),gamma=1.4,data=iClim,method="GCV.Cp")
      iClim$Tmid <- fit$residuals; #plot(fit,main="Tmid");
    }
    
    #make monthly averages 
    tempAggWT<-aggLags(datC=iClim, fun=mean, meas=Tmid, seg=month)
    subAgg<-tempAggWT[[1]]
    subAgg$month<-as.numeric(as.character(subAgg$month))
    subAgg<-subAgg[with(subAgg,order(year,month)),]
    
    tempWaveWT<-appLags(agDatC=subAgg, nUnits=nUnits, name=paste("month","Tmid", sep=""), fun=tempAggWT[[2]])
    tempWaveWT<-na.omit(tempWaveWT)
    
    monthStart=monthStart
    subsWT<-subset(tempWaveWT,tempWaveWT$month==monthStart) #choose start month
    subsWT<-subsWT[,-which(names(subsWT)=="month")] #Get rid of "month" column
    subsWT$year<-as.numeric(as.character(subsWT$year)) #make "year" readable to merge
    
    #reshape to show temp over the last nUnits
    labs<-c("year","t.00", paste("t.0",1:9, sep=""),
            paste("t.",10:nUnits, sep=""))
    names(subsWT) <- labs
    
    tempAggWP<-aggLags(datC=iClim, fun=mean, meas=Ppt, seg=month)
    subAgg<-tempAggWP[[1]]
    subAgg$month<-as.numeric(as.character(subAgg$month))
    subAgg<-subAgg[with(subAgg,order(year,month)),]
    
    tempWaveWP<-appLags(agDatC=subAgg, nUnits=nUnits, name=paste("month","Ppt", sep=""), fun=tempAggWP[[2]])
    tempWaveWP<-na.omit(tempWaveWP)
    subsWP<-subset(tempWaveWP,tempWaveWP$month==monthStart) #choose a start month
    subsWP<-subsWP[,-which(names(subsWP)=="month")] #remove month column
    subsWP$year<-as.numeric(as.character(subsWP$year)) #make "year" readable to merge
    
    labs<-c("year","p.00",paste("p.0",1:9, sep=""),
            paste("p.",10:nUnits, sep=""))
    names(subsWP) <- labs
    
    climDat<-cbind(subsWP,subsWT)
    
  }
  if(state=="Arizona"){ #quad-by-quad climate
    
    #Load Climate Data
    climDfile=paste(dataDir,"/climateData/monthClim.csv",sep="")
    iClim=read.csv(file=climDfile)
    iClim<-iClim[,-1]
    iClim$yearQ<-paste(iClim$year,iClim$quad, sep="_")
    
    #Detrend by month
    if(DETREND) {
      require(mgcv); #par(mfrow=c(2,2))
      
      fit=gam(Ppt~s(month,bs="cc")+factor(quad),gamma=1.4,data=iClim,method="GCV.Cp")
      iClim$Ppt <- fit$residuals; #plot(fit,main="Precip");
      
      fit=gam(Tmid~s(month,bs="cc")+factor(quad),gamma=1.4,data=iClim,method="GCV.Cp")
      iClim$Tmid <- fit$residuals; #plot(fit,main="Tmid");
    }
    
    #make monthly averages 
    tempAggWT<-aggLagsQ(datC=iClim, fun=mean, meas=Tmid, seg=month)
    subAgg<-tempAggWT[[1]]
    subAgg$month<-as.numeric(as.character(subAgg$month))
    subAgg<-subAgg[with(subAgg,order(yearQ,month)),]
    
    tempWaveWT<-appLags(agDatC=subAgg, nUnits=nUnits, name=paste("month","Tmid", sep=""), fun=tempAggWT[[2]])
    tempWaveWT<-na.omit(tempWaveWT)
    
    monthStart=monthStart
    subsWT<-subset(tempWaveWT,tempWaveWT$month==monthStart) #choose start month
    subsWT<-subsWT[,-which(names(subsWT)=="month")] #Get rid of "month" column
    elems<-unlist(strsplit(as.character(subsWT$yearQ) , "\\_" ) )
    m <- as.data.frame(matrix(elems , ncol = 2 , byrow = TRUE ))
    names(m)<-c("year","quad")
    subsWT<-cbind(m,subsWT[,-which(names(subsWT)=="yearQ")])
    subsWT$year<-as.numeric(as.character(subsWT$year)) #make "year" readable to merge
    
    #reshape to show temp over the last nUnits
    labs<-c("year","quad","t.00", paste("t.0",1:9, sep=""),
            paste("t.",10:nUnits, sep=""))
    names(subsWT) <- labs
    
    tempAggWP<-aggLagsQ(datC=iClim, fun=mean, meas=Ppt, seg=month)
    subAgg<-tempAggWP[[1]]
    subAgg$month<-as.numeric(as.character(subAgg$month))
    subAgg<-subAgg[with(subAgg,order(yearQ,month)),]
    
    tempWaveWP<-appLags(agDatC=subAgg, nUnits=nUnits, name=paste("month","Ppt", sep=""), fun=tempAggWP[[2]])
    tempWaveWP<-na.omit(tempWaveWP)
    subsWP<-subset(tempWaveWP,tempWaveWP$month==monthStart) #choose a start month
    subsWP<-subsWP[,-which(names(subsWP)=="month")] #remove month column
    elems<-unlist(strsplit(as.character(subsWP$yearQ) , "\\_" ) )
    m <- as.data.frame(matrix(elems , ncol = 2 , byrow = TRUE ))
    names(m)<-c("year","quad")
    subsWP<-cbind(m,subsWP[,-which(names(subsWP)=="yearQ")])
    subsWP$year<-as.numeric(as.character(subsWP$year)) #make "year" readable to merge
    
    labs<-c("year","quad","p.00",paste("p.0",1:9, sep=""),
            paste("p.",10:nUnits, sep=""))
    names(subsWP) <- labs
    
    climDat<-cbind(subsWP,subsWT)
    
  }
  
  
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
    
    ################# drop seedlings 
    #e <- (D$logarea.t0>0)
    #D <- D[e,];
    
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
    
    #colnames(sppData)<-c("logarea.t1","quad","year","age","logarea.t0","Group", colnames(sppData)[(1+nonCompLength.g):length(colnames(sppData))])
    
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
    
    #Drop seedlings
    #e <- (D$logarea>0)
    #D <- D[e,];
    
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

    ### Change this so sppData includes the response 
    #colnames(sppData)<-c("survives","age","ID","trackID","year","logarea","Group","quad",colnames(sppData)[(1+nonCompLength.s):length(colnames(sppData))])
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
    #year=as.numeric(as.factor(D$year))
    D$fyear = factor(D$year); 
    #Nyrs=length(unique(D$year))
    #Group=factor(D$Group)
    #Ngroups=length(unique(Group))
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
      dataG<-data.frame(intraD[,c("logarea.t1","year","ID","age","logarea.t0","Group","quad")], W)
      if(state=="Arizona"){climDat<-climDat[,-c(grep("year",names(climDat))[2],grep("quad",names(climDat))[2])]}else{
        climDat<-climDat[,-c(grep("year",names(climDat))[2])]
      }
      dataG$year<-as.numeric(paste(19,dataG$year,sep=""))
      
      if(state=="Arizona"){datag<-merge(dataG, climDat, by=c("year","quad"))}else{
        datag<-merge(dataG, climDat, by="year")}
      
      #datag$year<-as.factor(datag$year);
      datag$Group <- as.factor(datag$Group);
      
      ## define and enter covariates the way the fitting function wants
      pvars <- which(substr(names(datag),1,2)=="p."); 
      datag$pcovar <- as.matrix(datag[,pvars])
      
      tvars <- which(substr(names(datag),1,2)=="t."); 
      datag$tcovar <- as.matrix(datag[,tvars]) 
      
      lags <- matrix(0,nrow(datag),length(tvars)); 
      for(i in 1:ncol(lags)) lags[,i]=i; 
      datag$lags=as.matrix(lags); 
      
      data<-datag
      
    }
    if(surv==TRUE){
      
      dataS<-data.frame(intraD[,c("survives","year","ID","logarea","age","Group","quad")], W)
      
      if(state=="Arizona"){climDat<-climDat[,-c(grep("year",names(climDat))[2],grep("quad",names(climDat))[2])]}else{
        climDat<-climDat[,-c(grep("year",names(climDat))[2])]
      }
      dataS$year<-as.numeric(paste(19,dataS$year,sep=""))
      
      if(state=="Arizona"){datas<-merge(dataS, climDat, by=c("year","quad"))}else{
        datas<-merge(dataS, climDat, by="year")}
      #datas$year<-as.factor(datas$year);
      datas$Group <- as.factor(datas$Group);
      
      ## define and enter covariates the way the fitting function wants
      pvars <- which(substr(names(datas),1,2)=="p."); 
      datas$pcovar <- as.matrix(datas[,pvars])
      
      tvars <- which(substr(names(datas),1,2)=="t."); 
      datas$tcovar <- as.matrix(datas[,tvars]) 
      
      lags <- matrix(0,nrow(datas),length(tvars)); 
      for(i in 1:ncol(lags)) lags[,i]=i; 
      datas$lags=as.matrix(lags); 
      
      data<-datas
    }#end survival T/F
  }# W import for surv and grow T/F merge
  
  if(recr==TRUE){ #for recruitment, compute true parents
    dataR<-data.frame(intraD)
    
    if(state=="Arizona"){climDat<-climDat[,-c(grep("year",names(climDat))[2],grep("quad",names(climDat))[2])]}else{
      climDat<-climDat[,-c(grep("year",names(climDat))[2])]
    }
    
    dataR$year<-as.numeric(paste(19,dataR$year,sep=""))
    
    if(state=="Arizona"){datar<-merge(dataR, climDat, by=c("year","quad"))}else{
      datar<-merge(dataR, climDat, by="year")}
    
    #datas$year<-as.factor(datas$year);
    datar$Group <- as.factor(datar$Group);
    
    ## define and enter covariates the way the fitting function wants
    pvars <- which(substr(names(datar),1,2)=="p."); 
    datar$pcovar <- as.matrix(datar[,pvars])
    
    tvars <- which(substr(names(datar),1,2)=="t."); 
    datar$tcovar <- as.matrix(datar[,tvars]) 
    
    lags <- matrix(0,nrow(datar),length(tvars)); 
    for(i in 1:ncol(lags)) {lags[,i]=i;} 
    datar$lags=as.matrix(lags); 
    
    data<-datar
  } # and merge
  
  pvars <- which(substr(names(climDat),1,2)=="p."); 
  climDat$pcovar <- as.matrix(climDat[,pvars])
  
  tvars <- which(substr(names(climDat),1,2)=="t."); 
  climDat$tcovar <- as.matrix(climDat[,tvars]) 
  
  chuDat<-read.csv(paste(dataDir,"/climateData/Climate.csv",sep=""), header=T)
  
  return(list(data=data,climDat=climDat, chuDat=chuDat, iClim=iClim))
}


  
  