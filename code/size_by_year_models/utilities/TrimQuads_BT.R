##########################################################################
# This function removes individuals who would be flagged as seedlings
# based on the rule: (age==1)&(size <= log(.25)) ==> seedling 
# but who might not be seedlings because their age is > 1
#
# Individuals larger than 0.25 are assumed to be older (true age > 1) 
# regardless of what age is recorded for them. 
##########################################################################

trimQuadrats <- function(data,skip.years=TRUE,df,trans) {  
  
  if(trans=="surv"){data$logarea.t0=data$logarea}
  
  quadrats=unique(data$quad); doubtful=0; 
  for(j in 1:length(quadrats)) {
    qj = quadrats[j]; years = as.numeric(as.character(unique(data$year[data$quad == qj]))); 
    year1 = min(years)
    novice = which((data$year==year1)&(data$age==1)&(data$quad==qj)&(data$logarea.t0<log(0.25001)));  
    if(length(novice)>0) {
      data=data[-novice,];  #cat(qj,year1, length(novice), "\n"); 
      doubtful=doubtful+length(novice)
    } 
    
    
    if(skip.years){ # trim dubious ages after a gap in the data for a quadrat
      skipyr = which(diff(years)>1); 
      if(length(skipyr)>0) {
        skipyr=years[1+skipyr]; 
        for(k in 1:length(skipyr)){
          novice = which((data$year==skipyr[k])&(data$age==1)&(data$seedling==0)&(data$quad==qj)&(data$logarea.t0<log(0.25001)))  
          if(length(novice)>0) {
            data=data[-novice,];  #cat(qj,skipyr[k], length(novice), "\n"); 
            doubtful=doubtful+length(novice);
          }
        } 
      }
    }
  }
  likely.seedling=which((data$age==1)&(data$logarea.t0<=log(0.25001)))    
  data$year<-as.factor(data$year)
  
  # Drop likely seedlings or not 
  if(df=="df23"){
    if(length(likely.seedling)>0) data=data[-(likely.seedling), ]; 
  }
  
  if(df=="df1"){
    if(length(likely.seedling)>0) data=data[likely.seedling, ]; 
  }
  
  return(list(data=data)); 
}


