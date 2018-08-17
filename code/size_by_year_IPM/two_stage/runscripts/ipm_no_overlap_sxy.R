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
args       <- commandArgs(trailingOnly = F)
myargument <- args[length(args)]
myargument <- sub("-","",myargument)
i          <- as.numeric(myargument)



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
work_dir <- paste0(root,"/drivers/empirical/size_by_year_IPM/two_stage/")
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
grow_params_adults <- readRDS(file = grow_file)[[size_by_year]]
grow_params_kids   <- readRDS(file = grow_file)[["juveniles"]]

##  Load survival regression parameters
surv_file <- paste0("./vital_rate_params/surv_params_",state,"_",do_species,".RDS")
surv_params_adults <- readRDS(file = surv_file)[[size_by_year]]
surv_params_kids   <- readRDS(file = surv_file)[["juveniles"]]

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

spec_list  <- get_ipm_specs(grow_data)
inits_list <- get_initial_conditions(tlimit = spec_list[["tlimit"]],
                                     u = spec_list[["u"]],
                                     h = spec_list[["h"]],
                                     A = spec_list[["A"]],
                                     surv_data = surv_data)

##  Unlist IPM specifications and initial conditions to named objects
for(i in 1:length(names(spec_list))){
  assign(names(spec_list)[[i]], spec_list[[i]])
}
for(i in 1:length(names(inits_list))){
  assign(names(inits_list)[[i]], inits_list[[i]])
}

####
#### FIGURE OUT WHERE NEW RECRUITS GO IN THE nt1 VECTOR 
####
kidsize <- log(0.25)
kidbin <- 1 + sum(u<kidsize) 

if(transient == TRUE){
  start_kids <- nt1
  nt2[]      <- 0
}



####
####  ITERATE IPM ----
####
if(transient == TRUE){ tlimit <- 5000 }

##  Sample years
the_years <- Reduce(intersect, list(as.numeric(rec_params$year), 
                                    as.numeric(surv_params_adults$year_id), 
                                    as.numeric(grow_params_adults$year_id)))
all_years <- as.numeric(the_years)
year_save <- sample(all_years, tlimit, replace = TRUE)

##  Create empty data frame with necessary dimensions to store simulation
ipm_sims <- as.data.frame(matrix(NA, nrow = tlimit, ncol = 6)) %>%
  rename(iteration  = V1,
         year       = V2,
         cov_kids   = V3,
         cov_adults = V4,
         N_kids     = V5,
         N_adults   = V6) %>%
  mutate(iteration = 1:tlimit,
         year = year_save)

##  Plug in inital values
ipm_sims$cov_kids[1]   <- cov_save_kids[1]
ipm_sims$cov_adults[1] <- cov_save_adults[1]
ipm_sims$N_kids[1]     <- size_save_kids[1]
ipm_sims$N_adults[1]   <- size_save_adults[1]

size_save <- matrix(NA,length(u),tlimit)
tnt  <- nt1+nt2
size_save[,1] <- tnt/sum(tnt)

##  Calculate crowding variables
Wvars <- get_crowd_vars(data_dir = data_dir, 
                        state    = state, 
                        species  = do_species, 
                        b.r      = b.r)

##  Calculate the seedling growth size distribution; this is constant
grow_dist  <- get_growth_kids(u,u,grow_params_kids,kidsize)
if(round(sum(h*grow_dist))!=1){
  stop("Check grow_dist function, h*grow_dist should equal 1","\n")
}

##  Run IPM from t = 2 to t = tlimit, referenced by 'iiter'
pb <- txtProgressBar(min=2, max=tlimit, char="+", style=3, width=65)

for(iiter in 2:tlimit){
  ##  Set t0 and t1 for indexing; grab the first year
  t0      <- iiter - 1
  t1      <- iiter
  do_year <- ipm_sims$year[t1] # for random year effects
  
  ##  Calculate total basal cover for recs_per_area and crowding
  cover <- sum(ipm_sims$cov_kids[t0], ipm_sims$cov_adults[t0])

  ##  Calculate total population size distribution
  tnt  <- nt1+nt2
  expu <- exp(u)
  
  ##  Calculate crowding vector across sizes
  if(do_species %in% c("BOGR","SCSC")){
    Cr <- approxfun(b.r,h*c(0,cumsum(expu*tnt)),rule=2) # avoids integration issues
  }else{
    Cr <- splinefun(b.r,h*c(0,cumsum(expu*tnt)),method="natural")
  }
  
  Wvec <- makeW(Wvars,u,A,Cr)

  if(cover > 0){ # Don't do anything if species is extinct, for speed
    ##  Calculate recruits
    rpa      <- get_rec_number(params = rec_params, cover = cover, do_year = do_year)
    recruits <- f(v=u,u=u,Rpars=rec_size_params,rpa = rpa) 
    
    ##  Calculate growth*survival kernel / projection matrix
    # Kids/seedlings
    # growers_nt <- grow_dist * sumN(nt1,h) * get_survival_kids(u,u,Wvec,surv_params_adults,do_year)
    kid_survival_vec <- get_survival_kids(u,u,Wvec,surv_params_adults,do_year)
    growers_nt <- grow_dist * sumN(nt1,h) * kid_survival_vec[kidbin]
    
    # Adult plants
    Pmatrix_adults <- make.P.matrix(v=u, Wvec=Wvec, Gpars=grow_params_adults, Spars = surv_params_adults, do_year = do_year, G=G, S=S)
    
    ##  Project populations
    new_nt2         <- growers_nt + Pmatrix_adults%*%nt2
    new_nt1         <- rep(0,length(u)) 
    new_nt1[kidbin] <- recruits # put new recruits into their size bin, such that h*sum(new_nt1) = #recruits this year
    size_save[,iiter] <- (new_nt1+new_nt2)/sum((new_nt1+new_nt2))
  } # end if cover > 0
  
  ##  Update nt's for next timestep
  nt1 <- new_nt1 # kids
  nt2 <- new_nt2 # adults
  
  ##  Save population states
  ipm_sims$cov_kids[t1]   <- sumCover(u,nt1,h,A) # store the % cover as cm^2/cm^2
  ipm_sims$cov_adults[t1] <- sumCover(u,nt2,h,A) # store the % cover as cm^2/cm^2
  ipm_sims$N_kids[t1]     <- sumN(nt1,h) # sum abundance over size stages
  ipm_sims$N_adults[t1]   <- sumN(nt2,h) # sum abundance over size stages
  
  if(transient == TRUE){
    if(iiter %in% seq(1,tlimit,100)) {
      nt1 <- start_kids # kids
      nt2[] <- 0 # adults
      
      ##  Save population states
      ipm_sims$cov_kids[t1]   <- sumCover(u,nt1,h,A) # store the % cover as cm^2/cm^2
      ipm_sims$cov_adults[t1] <- sumCover(u,nt2,h,A) # store the % cover as cm^2/cm^2
      ipm_sims$N_kids[t1]     <- sumN(nt1,h) # sum abundance over size stages
      ipm_sims$N_adults[t1]   <- sumN(nt2,h) # sum abundance over size stages
    }
  }
  
  setTxtProgressBar(pb, iiter)
} # end IPM iterations



####
####  SAVE OUTPUT AND PLOT (OPTIONAL) ----
####
out_file <- paste0(out_dir, paste(state, do_species, size_by_year, sep = "_"),".RDS")
saveRDS(ipm_sims, out_file)

if(size_by_year == "yearonly"){
  dir.create(file.path("../ipm_results/size_distributions/"), showWarnings = FALSE) # create a folder for results
  saveRDS(size_save[,(burn_in+1):ncol(size_save)], paste0("../ipm_results/size_distributions/",paste(state, do_species, "sizeDist", sep = "_"),".RDS"))
}
