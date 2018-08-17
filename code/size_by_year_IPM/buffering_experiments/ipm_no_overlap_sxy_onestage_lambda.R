################################################################################
##  ipm_no_overlap_sxy_onstage_lambda.R: This is the typical one-stage IPM, set
##  up to calculate the low density growth rate of each species under different
##  size*year perturbations. The script is set up to be run by the bash script
##  'run_lambda_experiments.sh' in this same directory.
##
##  ____________________________________________________________________________
##  Author: Andrew Tredennick (atredenn@gmail.com)
##  Contributors: Peter Adler, Britta Teller, Stephen Ellner
##  Date created: November 28, 2017
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
isim <- as.numeric(myargument)



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
out_dir  <- paste0(root,"/drivers/empirical/size_by_year_IPM/buffering_experiments/lambda_results/")

setwd(work_dir) # set working directory
dir.create(file.path(out_dir), showWarnings = FALSE) # create a folder for results



####
####  SOURCE ALL NEEDED FUNCTIONS ----
####
files_to_source <- list.files("../buffering_experiments/ipm_onestage_functions/")
for(dofile in files_to_source){
  source(paste0("../buffering_experiments/ipm_onestage_functions/",dofile))
}
source("../one_stage/vital_rate_functions/fetch_data_fxns.R")
source("../buffering_experiments/buffering_fxns.R")
source("../one_stage/vital_rate_functions/reformat_state_df_fxns.R")



####
####  STATE, SPECIES, AND SIMULATION INFORMATION ----
####
# Make matrix of state-species for one-stage modeling
state_species <- data.frame(state = c(rep("montana",4),
                                      rep("arizona",1),
                                      rep("newmexico",1),
                                      rep("kansas",3)),
                            species = c(sort(c("BOGR","HECO","PASM","POSE")),
                                        sort(c("BORO")),
                                        sort(c("BOER")),
                                        sort(c("ANGE","BOCU","BOHI"))))



####
####  PERFORM ANALYSIS, LOOPING OVER SPECIES ----
####
for(doit in 1:nrow(state_species)){
  state      <- as.character(state_species[doit, "state"])
  do_species <- as.character(state_species[doit, "species"])
  
  ## Run it all...
  
  transform_matrix <- expand.grid(c("variable_small","variable_large","sign_flip"),
                                  c("variable_small","variable_large","sign_flip"))
  transform_matrix <- rbind(transform_matrix,
                            data.frame(Var1 = c("variable_small","variable_large","sign_flip"),
                                       Var2 = rep("none",3)),
                            data.frame(Var1 = rep("none",3),
                                       Var2 = c("variable_small","variable_large","sign_flip")),
                            data.frame(Var1 = "none", Var2 = "none"))
  colnames(transform_matrix) <- c("growth","survival")
  grow_transform <- as.character(transform_matrix[isim,"growth"])
  surv_transform <- as.character(transform_matrix[isim,"survival"])
  
  
  
  ####
  ####  GET STABLE SIZE DISTRIBUTION ----
  ####
  fname <- paste0("../ipm_results/size_distributions/",state,"_",do_species,"_sizeDist.RDS")
  size_dist <- readRDS(fname)
  stable_size <- rowMeans(size_dist)
  fix_cover <- 0.001
  
  
  
  ####
  ####  LOAD VITAL RATE PARAMETERS ----
  ####
  ##  Load growth regression parameters
  grow_file <- paste0("./vital_rate_params/growth_params_",state,"_",do_species,".RDS")
  grow_params <- readRDS(file = grow_file)[["yearonly"]]
  
  ##  Load survival regression parameters
  surv_file <- paste0("./vital_rate_params/surv_params_",state,"_",do_species,".RDS")
  surv_params <- readRDS(file = surv_file)[["yearonly"]]
  
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
  
  ##  Fit growth and survival models
  gdata <- fetch_grow_data(data_dir, state, do_species)
  gfit  <- lmer(logarea.t1 ~ logarea.t0 * W + (1|Group) + (1|fyear), 
                data = reformat_growth(gdata,state))
  sdata <- fetch_surv_data(data_dir, state, do_species)
  sfit  <- glmer(survives ~ logarea * W + (1|Group) + (1|fyear), 
                 data = reformat_survival(sdata,state), family = "binomial")
  
  ##  Unlist IPM specifications to named objects
  spec_list  <- get_ipm_specs(grow_data)
  for(i in 1:length(names(spec_list))){
    assign(names(spec_list)[[i]], spec_list[[i]])
  }
  
  
  
  ####
  ####  DEFINE SIZE-BY-YEAR EFFECTS ----
  ####
  ##  Growth
  if(grow_transform == "variable_small"){
    grow_size_year <- make_small_variable(sizes = u, mod_fit = gfit, view_plot = F)
  }
  if(grow_transform == "variable_large"){
    grow_size_year <- make_large_variable(sizes = u, mod_fit = gfit, view_plot = F)
  }
  if(grow_transform == "sign_flip"){
    grow_size_year <- make_sign_flip(sizes = u, mod_fit = gfit, view_plot = F)
  }
  if(grow_transform == "none"){
    years <- as.numeric(row.names(coef(gfit)$fyear))
    grow_size_year <- matrix(0,length(years),length(u))
    row.names(grow_size_year) <- years
  }
  
  ##  Survival
  if(surv_transform == "variable_small"){
    surv_size_year <- make_small_variable(sizes = u, mod_fit = sfit, view_plot = F)
  }
  if(surv_transform == "variable_large"){
    surv_size_year <- make_large_variable(sizes = u, mod_fit = sfit, view_plot = F)
  }
  if(surv_transform == "sign_flip"){
    surv_size_year <- make_sign_flip(sizes = u, mod_fit = sfit, view_plot = F)
  }
  if(surv_transform == "none"){
    years <- as.numeric(row.names(coef(sfit)$fyear))
    surv_size_year <- matrix(0,length(years),length(u))
    row.names(surv_size_year) <- years
  }
  
  
  
  # set to stable size distribution
  nt <- stable_size
  # set to fix cover value
  tmp <- fix_cover*100/(h*sum(nt*exp(u)))
  nt <- nt*tmp
  
  
  
  ####
  ####  ITERATE IPM ----
  ####
  ##  Create year vector for observed years
  the_years <- Reduce(intersect, list(rec_params$year, surv_params$year_id, grow_params$year_id))
  all_years <- as.numeric(the_years)
  year_save <- sort(all_years)
  tlimit <- length(year_save)+1
  
  ##  Create empty data frame with necessary dimensions to store simulation
  ipm_sims <- as.data.frame(matrix(NA, nrow = tlimit, ncol = 4)) %>%
    rename(iteration = V1, year = V2, cover = V3, abundance = V4) %>%
    mutate(iteration = 1:tlimit, year = c(1,year_save))
  
  ##  Plug in initial values
  ipm_sims$cover[1]     <- sumCover(u,nt,h,A)
  ipm_sims$abundance[1] <- sumN(nt,h)
  
  ##  Calculate crowding variables
  Wvars <- get_crowd_vars(data_dir = data_dir, 
                          state    = state, 
                          species  = do_species, 
                          b.r      = b.r)
  
  max_size <- log(maxSize)
  fecundity_cutoff <- sum(u<max_size) 
  
  for(iiter in 2:tlimit){
    ##  Set t0 and t1 for indexing; grab the first year
    t0      <- iiter - 1
    t1      <- iiter
    do_year <- ipm_sims$year[t1] # for random year effects
    
    ##  Calculate total basal cover for recs_per_area and crowding
    cover <- fix_cover/100
    
    ##  Calculate total population size
    N    <- sumN(nt,h)
    expu <- exp(u)
    
    ##  Calculate crowding vector across sizes
    if(do_species %in% c("BOGR","BORO")){
      Cr <- approxfun(b.r,h*c(0,cumsum(expu*nt)),rule=2) # avoids integration issues
    }else{
      Cr <- splinefun(b.r,h*c(0,cumsum(expu*nt)),method="natural")
    }
    
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
                               surv_allsize_years = surv_size_year,
                               grow_allsize_years = grow_size_year,
                               fecundity_cutoff = fecundity_cutoff)	
      
      ##  Project population
      new_nt <- Kmatrix%*%nt
      
    } # end if cover > 0
    
    ##  Update nt's for next timestep
    nt <- new_nt
    
    ##  Save population states
    ipm_sims$cover[t1] <- sumCover(u,nt,h,A)
    ipm_sims$abundance[t1] <- sumN(nt,h)
    
    ##  Reset to low cover
    tmp <- fix_cover*100/(h*sum(nt*exp(u)))
    nt  <- nt*tmp
  } # end IPM iterations
  
  
  
  ####
  ####  SAVE OUTPUT ----
  ####
  out_file <- paste0(out_dir, paste(state, do_species, isim, sep = "_"),".RDS")
  saveRDS(ipm_sims, out_file)
} # End state-species loop



