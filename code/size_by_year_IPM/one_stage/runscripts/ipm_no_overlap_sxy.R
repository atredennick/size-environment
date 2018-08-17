################################################################################
##  ipm_no_overlap_sxy.R: R script that can be called from bash script
##  "run_ipms.sh" to run single species Integral Projection Model (IPM) 
##  simulations for each of our 15 focal species. IPMs are run with and without 
##  size-by-year interactions to analyze the effect of size-specific climate 
##  responses on long-term population dynamics. The code returns and saves a 
##  data frame of species cover over time, for each species and each model type 
##  (interaction or not). This code requires several utility functions from the 
##  './ipm_functions/' folder in this directory.
##
##  ____________________________________________________________________________
##  Author: Andrew Tredennick (atredenn@gmail.com)
##  Contributors: Peter Adler, Britta Teller, Stephen Ellner
##  Date created: August 28, 2017
################################################################################


##  Clear everything
rm(list = ls(all.names = TRUE))

####
####  GET SYSTEM ARGUMENTS ----
####
##  This is for running in batch later...set manually below for now
args <- commandArgs(trailingOnly = F)
myargument <- args[length(args)]
myargument <- sub("-","",myargument)
i <- as.numeric(myargument)

##  Set model specification manually using 'i'
##  i = 1: idaho/ARTR/size by year interaction
##  i = 30: kansas/bohi/size only
## i <- 1 # hard set for testing



####
####  LOAD LIBRARIES ----
####
library(tidyverse)
library(dplyr)
library(gamlss)
library(stats4)
library(boot) 
library(lme4)


####
####  SET DIRECTORIES ----
####
root     <- ifelse(.Platform$OS.type=="windows","c:/repos","~/repos")
work_dir <- paste0(root,"/drivers/empirical/size_by_year_IPM/one_stage/")
data_dir <- paste0(root,"/driversdata/data/")
out_dir  <- paste0(root,"/drivers/empirical/size_by_year_IPM/ipm_results/equilibrium_runs/")

setwd(work_dir) # set working directory
dir.create(file.path(out_dir), showWarnings = FALSE) # create a folder for results



####
####  SOURCE ALL NEEDED FUNCTIONS ----
####
files_to_source <- list.files("./ipm_functions/")
for(dofile in files_to_source){
  source(paste0("./ipm_functions/",dofile))
}



####
####  STATE, SPECIES, AND SIMULATION INFORMATION ----
####
state_species <- data.frame(state = c(rep("idaho",4),
                                      rep("montana",4),
                                      rep("arizona",2),
                                      rep("newmexico",2),
                                      rep("kansas",3)),
                            species = c(sort(c("ARTR","HECO","POSE","PSSP")),
                                        sort(c("BOGR","HECO","PASM","POSE")),
                                        sort(c("BOER","BORO")),
                                        sort(c("BOER","SPFL")),
                                        sort(c("ANGE","BOCU","BOHI"))))
state_species$size_by_year <- "sizeXyear"
state_species <- rbind(state_species,
                       mutate(state_species, size_by_year="yearonly"))
state_species$transient <- FALSE
state_species <- rbind(state_species,
                       mutate(state_species, transient=TRUE))

state        <- as.character(state_species[i,"state"])
do_species   <- as.character(state_species[i,"species"])
size_by_year <- state_species[i,"size_by_year"]
transient    <- state_species[i,"transient"]

if(transient == TRUE){
  out_dir  <- paste0(root,"/drivers/empirical/size_by_year_IPM/ipm_results/transient_runs/")
  dir.create(file.path(out_dir), showWarnings = FALSE) # create a folder for results
}



####
####  LOAD VITAL RATE PARAMETERS ----
####
##  Load growth regression parameters
grow_file <- paste0("./vital_rate_params/growth_params_",state,"_",do_species,".RDS")
grow_params <- readRDS(file = grow_file)[[size_by_year]]

##  Load survival regression parameters
surv_file <- paste0("./vital_rate_params/surv_params_",state,"_",do_species,".RDS")
surv_params <- readRDS(file = surv_file)[[size_by_year]]

##  Load recruitment regression parameters
rec_file <- paste0("./vital_rate_params/rec_params_",state,"_",do_species,".RDS")
rec_params <- readRDS(file = rec_file)
rec_size_params <- get_rec_size(data_dir,state,do_species)



####
####  SET UP IPM SPECIFICATIONS AND INITIAL CONDITIONS ----
####
if(state %in% c("arizona", "kansas", "newmexico")){
  grow_data <- read.csv(paste0(data_dir,state,"/speciesData/",do_species,"/growDnoNA_reduced.csv"))
  surv_data <- read.csv(paste0(data_dir,state,"/speciesData/",do_species,"/survD_reduced.csv"))
}
if(state %in% c("idaho", "montana")){
  grow_data <- read.csv(paste0(data_dir,state,"/speciesData/",do_species,"/growDnoNA.csv"))
  surv_data <- read.csv(paste0(data_dir,state,"/speciesData/",do_species,"/survD.csv"))
}

##  Unlist IPM specifications to named objects
spec_list  <- get_ipm_specs(grow_data)
for(i in 1:length(names(spec_list))){
  assign(names(spec_list)[[i]], spec_list[[i]])
}



####
#### FIGURE OUT WHERE NEW RECRUITS GO IN THE nt1 VECTOR 
####
if(transient == TRUE){
  kidsize <- log(0.25)
  kidbin <- 1 + sum(u<kidsize) 
  nt[] <- 0
  nt[kidbin] <- 94
  start_kids <- nt
}


####
####  ITERATE IPM ----
####
if(transient == TRUE){ tlimit <- 5000 } 

##  Sample years
the_years <- Reduce(intersect, list(rec_params$year, surv_params$year_id, grow_params$year_id))
all_years <- as.numeric(the_years)
year_save <- sample(all_years, tlimit, replace = TRUE)

##  Create empty data frame with necessary dimensions to store simulation
ipm_sims <- as.data.frame(matrix(NA, nrow = tlimit, ncol = 4)) %>%
  rename(iteration = V1, year = V2, cover = V3, abundance = V4) %>%
  mutate(iteration = 1:tlimit, year = year_save)

##  Plug in initial values
ipm_sims$cover[1]     <- sumCover(u,nt,h,A)
ipm_sims$abundance[1] <- sumN(nt,h)

size_save <- matrix(NA,length(u),tlimit)
size_save[,1] <- nt/sum(nt)

##  Calculate crowding variables
Wvars <- get_crowd_vars(data_dir = data_dir, 
                        state    = state, 
                        species  = do_species, 
                        b.r      = b.r)

max_size <- log(maxSize)
fecundity_cutoff <- sum(u<max_size) 

##  Run IPM from t = 2 to t = tlimit, referenced by 'iiter'
pb <- txtProgressBar(min=2, max=tlimit, char="+", style=3, width=65) # progress bar

for(iiter in 2:tlimit){
  
  ##  Set t0 and t1 for indexing; grab the first year
  t0      <- iiter - 1
  t1      <- iiter
  do_year <- ipm_sims$year[t1] # for random year effects
  
  ##  Calculate total basal cover for recs_per_area and crowding
  cover <- ipm_sims$cover[t0]

  ##  Calculate total population size
  N    <- sumN(nt,h)
  expu <- exp(u)

  ##  Calculate crowding vector across sizes
  Cr <- splinefun(b.r,h*c(0,cumsum(expu*nt)),method="natural")
  # Cr <- approxfun(b.r,h*c(0,cumsum(expu*nt)),rule=2) # avoids integration issues
  
  ##  Calculate crowding vector
  Wvec <- makeW(Wvars,u,A,Cr)
  
  ##  Calculate recruits per area
  rpa <- get_recs_per_area(params = rec_params, cover = cover, A = A, do_year = do_year)

  if(cover > 0){ # Don't do anything if species is extinct, for speed
    
    Kmatrix <- make.K.matrix(v       = u,
                             muW     = Wvec,
                             Rpars   = rec_size_params,
                             rpa     = rpa,
                             Gpars   = grow_params,
                             Spars   = surv_params,
                             do_year = do_year,
                             fecundity_cutoff = fecundity_cutoff)	
    
    ##  Project population
    new_nt <- Kmatrix%*%nt
    size_save[,iiter] <- new_nt/sum(new_nt)

  } # end if cover > 0
  
  ##  Update nt's for next timestep
  nt <- new_nt
  
  ##  Save population states
  ipm_sims$cover[t1] <- sumCover(u,nt,h,A)
  ipm_sims$abundance[t1] <- sumN(nt,h)
  
  if(transient == TRUE){
    if(iiter %in% seq(1,tlimit,100)) {
      nt <- start_kids # kids
      
      ##  Save population states
      ipm_sims$cover[t1] <- sumCover(u,nt,h,A)
      ipm_sims$abundance[t1] <- sumN(nt,h)
    }
  }
  
  ##  Update progess
  setTxtProgressBar(pb, iiter)
} # end IPM iterations



####
####  SAVE OUTPUT ----
####
out_file <- paste0(out_dir, paste(state, do_species, size_by_year, sep = "_"),".RDS")
saveRDS(ipm_sims, out_file)

if(size_by_year == "yearonly"){
  dir.create(file.path("../ipm_results/size_distributions/"), showWarnings = FALSE) # create a folder for results
  saveRDS(size_save[,(burn_in+1):ncol(size_save)], paste0("../ipm_results/size_distributions/",paste(state, do_species, "sizeDist", sep = "_"),".RDS"))
}

