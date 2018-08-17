################################################################################
##  reformat_state_df_fxns.R: functions for reformatting dataframes with correct
##  columns for each state.
##
##  ____________________________________________________________________________
##  Author: Andrew Tredennick
##  Date created: 8-29-2017
################################################################################

#' Reformat growth data frames.
#' @author Andrew Tredennick
#' @param df Data frame of growth observations.
#' @param state Character scalar indicating the state we are dealing with.
#' @return A reformatted dataframe with correct columns.
reformat_growth <- function(df, state){
  D <- df
  do_site <- state
  # Add group info for each site individually, as needed
  if(do_site=="idaho"){ D <- D}
  
  if(do_site=="arizona"){
    D$Group=as.factor(substr(D$quad,1,1))
    D=subset(D,as.numeric(as.character(year))>=17)
  }
  
  if(do_site=="kansas")
    D$Group=as.numeric(D$Group)-1
  
  if(do_site=="montana"){
    # Remove suspect quadrat years
    # quadyears <- with(D, paste0(quad, year))
    # quadyrs_to_remove <- read.csv("../data/Montana/BOGR_edited/suspect_BOGR_quads.csv")
    # quadyrs_to_remove$year <- quadyrs_to_remove$year - 1900
    # bad_quadyears <- with(quadyrs_to_remove, paste0(quadrat,year))
    # torms <- which(quadyears %in% bad_quadyears)
    # if(length(torms)>0) { D <- D[-torms,] }
    
    # Then we moved some specific points:
    tmp2<-which(D$quad=="A12" & D$year==44)
    tmp3<-which(D$quad=="B1"  & D$year==44)
    tmp41<-which(D$quad=="E4" & D$year==33) 
    tmp42<-which(D$quad=="E4" & D$year==34) 
    tmp43<-which(D$quad=="E4" & D$year==43)
    tmp44<-which(D$quad=="E4" & D$year==44)
    tmpONE<-c(tmp2,tmp3,tmp41,tmp42,tmp43,tmp44)
    if(length(tmpONE)>0) D<-D[-tmpONE,]
    D$Group=as.factor(substr(D$quad,1,1)) 
  }
  
  if(do_site=="newmexico")
    D$Group=as.factor(substr(D$quad,1,1))
  
  # Get the years right for Kansas
  if(do_site=="kansas"){
    D$year <- as.numeric(as.character(D$year))
    D <- subset(D, year<68)
    ##to remove some points:
    #for q25
    tmp1<-which(D$quad=="q25")
    #for q27
    tmp2<-which(D$quad=="q27")
    #for q28
    tmp3<-which(D$quad=="q28")
    #for q30
    tmp4<-which(D$quad=="q30")
    #for q31
    tmp5<-which(D$quad=="q31" & (D$year<35 | D$year>39))
    #for q32
    tmp6<-which(D$quad=="q32" & (D$year<35 | D$year>41))
    tmp<-c(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)
    D<-D[-tmp,]
  }
  
  return(D)
}



#' Reformat survival data frames.
#' @author Andrew Tredennick
#' @param df Data frame of growth observations.
#' @param state Character scalar indicating the state we are dealing with.
#' @return A reformatted dataframe with correct columns.
reformat_survival <- function(df, state){
  D <- df
  do_site <- state
  # Add group info for each site individually, as needed
  if(do_site=="idaho"){ D <- D}
  
  # Add group info for each site individually, as needed
  if(do_site=="arizona"){
    D$Group=as.factor(substr(D$quad,1,1))
    D=subset(D,as.numeric(as.character(year))>=17)
  }
  if(do_site=="kansas")
    D$Group=as.numeric(D$Group)-1
  if(do_site=="montana"){
    # Remove suspect quadrat years
    # quadyears <- with(D, paste0(quad, year))
    # quadyrs_to_remove <- read.csv("../data/Montana/BOGR_edited/suspect_BOGR_quads.csv")
    # quadyrs_to_remove$year <- quadyrs_to_remove$year - 1900
    # bad_quadyears <- with(quadyrs_to_remove, paste0(quadrat,year))
    # torms <- which(quadyears %in% bad_quadyears)
    # if(length(torms)>0) { D <- D[-torms,] }
    
    # Then we moved some specific points:
    tmp2<-which(D$quad=="A12" & D$year==44)
    tmp3<-which(D$quad=="B1"  & D$year==44)
    tmp41<-which(D$quad=="E4" & D$year==33) 
    tmp42<-which(D$quad=="E4" & D$year==34) 
    tmp43<-which(D$quad=="E4" & D$year==43)
    tmp44<-which(D$quad=="E4" & D$year==44)
    tmpONE<-c(tmp2,tmp3,tmp41,tmp42,tmp43,tmp44)
    if(length(tmpONE)>0) D<-D[-tmpONE,]
    D$Group=as.factor(substr(D$quad,1,1)) 
  }
  if(do_site=="newmexico")
    D$Group=as.factor(substr(D$quad,1,1))
  
  # Get the years right for Kansas
  if(do_site=="kansas"){
    D$year <- as.numeric(as.character(D$year))
    D <- subset(D, year<68)
    ##to remove some points:
    #for q25
    tmp1<-which(D$quad=="q25")
    #for q27
    tmp2<-which(D$quad=="q27")
    #for q28
    tmp3<-which(D$quad=="q28")
    #for q30
    tmp4<-which(D$quad=="q30")
    #for q31
    tmp5<-which(D$quad=="q31" & (D$year<35 | D$year>39))
    #for q32
    tmp6<-which(D$quad=="q32" & (D$year<35 | D$year>41))
    tmp<-c(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)
    D<-D[-tmp,]
  }
  
  return(D)
}

