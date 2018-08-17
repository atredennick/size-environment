#' Fetches and formats growth data.
#' 
#' @author Andrew Tredennick
#' @param data_dir Path to directory containing data.
#' @param state  The state we are working with.
#' @param do_species The species we are working with.
#' @return A dataframe.
fetch_grow_data <- function(data_dir, state, do_species) {
  require(tidyverse)
  require(dplyr)
  
  ##  Read in demographic data
  if(state == "idaho" | state == "montana"){
    grow_file <- paste0(data_dir,state,"/speciesData/",do_species,"/growDnoNA.csv")
  }
  if(state == "arizona" | state == "newmexico" | state == "kansas"){
    grow_file <- paste0(data_dir,state,"/speciesData/",do_species,"/growDnoNA_reduced.csv")
  }
  
  ring_file <- paste0(data_dir,state,"/speciesdata/",do_species,"/",do_species,"_nbhood_rings.csv")
  
  ring_data <- read.csv(ring_file) %>%
    dplyr::rename(trackID = genetID) %>%
    dplyr::select(species, quad, year, trackID, starts_with(do_species))
  
  suppressWarnings( # safely ignore quad coerce to character warning
    growdata_full <- read.csv(grow_file) %>%
      dplyr::filter(allEdge==0) %>%
      dplyr::mutate(logarea.t0 = log(area.t0),
                    logarea.t1 = log(area.t1),
                    fyear      = factor(year)) %>%
      dplyr::left_join(ring_data, by = c("quad","year","trackID")) %>%
      dplyr::arrange(X)
  )
  
  ##  Create crowding effect
  distance_file <- paste0(data_dir,state,"/speciesData/",state,"DistanceWeights.csv")
  dist_wts <- read.csv(distance_file) %>%
    gather(key = Species, value = weight, -midRings) %>%
    filter(Species == do_species) %>%
    spread(key = Species, value = weight)
  
  comp_matrix <- data.matrix(dplyr::select(growdata_full, starts_with(do_species))) # matrix of conspecific areas in the annuli 
  wvec <- comp_matrix%*%dist_wts[,as.character(do_species)]
  
  growdata_full <- growdata_full %>%
    dplyr::mutate(W = as.numeric(wvec)) %>%
    dplyr::select(-starts_with(as.character(do_species)))
  
  if(do_species == "ARTR"){ # Remove outliers (large plants that obviously do 
                            # not turn into tiny plants) -- PBA
    growdata_full <- growdata_full %>%
      dplyr::filter(quad != "Q23" & year != 1945 & trackID != 67) %>%
      dplyr::filter(quad != "Q12" & year != 1955 & trackID != 25) %>%
      dplyr::filter(quad != "Q26" & year != 1945 & trackID != 73)
  }
  
  return(growdata_full)
}



#' Fetches and formats survival data.
#' 
#' @author Andrew Tredennick
#' @param data_dir Path to directory containing data.
#' @param state  The state we are working with.
#' @param do_species The species we are working with.
#' @return A dataframe
fetch_surv_data <- function(data_dir, state, do_species) {
  require(tidyverse)
  require(dplyr)
  
  ##  Read in demographic data
  if(state != "arizona" | state != "newmexico" | state != "kansas"){
    surv_file <- paste0(data_dir,state,"/speciesData/",do_species,"/survD.csv")
  }
  if(state == "arizona" | state == "newmexico" | state == "kansas"){
    surv_file <- paste0(data_dir,state,"/speciesData/",do_species,"/survD_reduced.csv")
  }
  ring_file <- paste0(data_dir,state,"/speciesdata/",do_species,"/",do_species,"_nbhood_rings.csv")
  
  ring_data <- read.csv(ring_file) %>%
    dplyr::rename(trackID = genetID) %>%
    dplyr::select(species, quad, year, trackID, starts_with(do_species))
  
  suppressWarnings( # safely ignore quad coerce to character warning
    survdata_full <- read.csv(surv_file) %>%
      dplyr::filter(allEdge == 0) %>%
      dplyr::mutate(logarea = log(area),
                    fyear   = factor(year)) %>%
      dplyr::left_join(ring_data, by = c("quad","year","trackID")) %>%
      dplyr::arrange(X)
  )
  
  ##  Create crowding effect
  distance_file <- paste0(data_dir,state,"/speciesData/",state,"DistanceWeights.csv")
  dist_wts <- read.csv(distance_file) %>%
    gather(key = Species, value = weight, -midRings) %>%
    filter(Species == do_species) %>%
    spread(key = Species, value = weight)
  
  comp_matrix <- data.matrix(dplyr::select(survdata_full, starts_with(do_species))) # matrix of conspecific areas in the annuli 
  wvec <- comp_matrix%*%dist_wts[,as.character(do_species)]
  
  survdata_full <- survdata_full %>%
    dplyr::mutate(W = wvec) %>%
    dplyr::select(-starts_with(as.character(do_species)))
  
  return(survdata_full)
}



#' Fetches and formats recruiment data
#' 
#' @author Andrew Tredennick
#' @param data_dir Path to directory containing data.
#' @param state  The state we are working with.
#' @param do_species The species we are working with.
#' @return A dataframe.
fetch_recruit_data <- function(data_dir, state, do_species) {
  require(tidyverse)
  require(dplyr)
  
  ##  Read in recruitment data
  if(state != "arizona" | state != "newmexico" | state != "kansas"){
    rec_file <- paste0(data_dir,state,"/speciesData/",do_species,"/recArea.csv")
  }
  if(state == "arizona" | state == "newmexico" | state == "kansas"){
    rec_file <- paste0(data_dir,state,"/speciesData/",do_species,"/recArea_reduced.csv")
  }
  
  tmpD <- read.csv(rec_file) 
  
  ##  Reformate data frame by state
  if(state=="arizona"){
    tmpD$Group=as.factor(substr(tmpD$quad,1,1))
    tmpD=subset(tmpD,year>=17)
  }
  if(state=="kansas")
    tmpD$Group=as.numeric(tmpD$Group)-1
  if(state=="montana"){
    # Remove suspect quadrat years
    # quadyears <- with(tmpD, paste0(quad, year))
    # quadyrs_to_remove <- read.csv("../data/Montana/BOGR_edited/suspect_BOGR_quads.csv")
    # quadyrs_to_remove$year <- quadyrs_to_remove$year - 1900
    # bad_quadyears <- with(quadyrs_to_remove, paste0(quadrat,year))
    # torms <- which(quadyears %in% bad_quadyears)
    # if(length(torms)>0) { tmpD <- tmpD[-torms,] }
    
    # Then we moved some specific points:
    tmp2<-which(tmpD$quad=="A12" & tmpD$year==44)
    tmp3<-which(tmpD$quad=="B1"  & tmpD$year==44)
    tmp41<-which(tmpD$quad=="E4" & tmpD$year==33) 
    tmp42<-which(tmpD$quad=="E4" & tmpD$year==34) 
    tmp43<-which(tmpD$quad=="E4" & tmpD$year==43)
    tmp44<-which(tmpD$quad=="E4" & tmpD$year==44)
    tmpONE<-c(tmp2,tmp3,tmp41,tmp42,tmp43,tmp44)
    if(length(tmpONE)>0) tmpD<-tmpD[-tmpONE,]
    tmpD$Group=as.factor(substr(tmpD$quad,1,1))
  }
  if(state=="newmexico")
    tmpD$Group=as.factor(substr(tmpD$quad,1,1))
  
  # Get the years right for Kansas
  if(state=="kansas"){
    tmpD <- subset(tmpD, year<68)
    ##to remove some points:
    #for q25
    tmp1<-which(tmpD$quad=="q25")
    #for q27
    tmp2<-which(tmpD$quad=="q27")
    #for q28
    tmp3<-which(tmpD$quad=="q28")
    #for q30
    tmp4<-which(tmpD$quad=="q30")
    #for q31
    tmp5<-which(tmpD$quad=="q31" & (tmpD$year<35 | tmpD$year>39))
    #for q32
    tmp6<-which(tmpD$quad=="q32" & (tmpD$year<35 | tmpD$year>41))
    tmp<-c(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)
    tmpD<-tmpD[-tmp,]
  }
  
  rec_data <- tmpD %>%
    dplyr::select(quad, year, NRquad, totParea, Group) %>%
    dplyr::mutate(species = as.character(do_species)) %>%
    dplyr::rename(recruits = NRquad, cover = totParea)
  
  ##  Calculate mean cover by group and year
  group_data <- rec_data %>%
    dplyr::group_by(year, Group) %>%
    dplyr::summarise(Gcov = mean(cover))
  
  ##  Merge the group level data back in
  rec_data <- rec_data %>%
    dplyr::left_join(group_data, by = c("year", "Group"))
  
  return(rec_data)
}

