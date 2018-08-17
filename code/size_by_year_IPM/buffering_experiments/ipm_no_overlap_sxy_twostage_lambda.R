################################################################################
##  ipm_no_overlap_sxy_twostage_lambda.R: This is the two-stage IPM, set up to
##  calculate the low density growth rate of each species under different
##  size*year perturbations. The script is set up to be run by the bash script
##  'run_lambda_experiments.sh' in this same directory.
##
##  ____________________________________________________________________________
##  Author: Andrew Tredennick (atredenn@gmail.com)
##  Contributors: Peter Adler, Britta Teller, Stephen Ellner
##  Date created: December 8, 2017
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
work_dir <- paste0(root,"/drivers/empirical/size_by_year_IPM/two_stage/")
data_dir <- paste0(root,"/driversdata/data/")
out_dir  <- paste0(root,"/drivers/empirical/size_by_year_IPM/buffering_experiments/lambda_results/")

setwd(work_dir) # set working directory
dir.create(file.path(out_dir), showWarnings = FALSE) # create a folder for results



####
####  SOURCE ALL NEEDED FUNCTIONS ----
####
files_to_source <- list.files("../buffering_experiments/ipm_twostage_functions/")
for(dofile in files_to_source){
  source(paste0("../buffering_experiments/ipm_twostage_functions/",dofile))
}
files_to_source <- list.files("../two_stage/vital_rate_functions/")
dont_source     <- grep("jags", files_to_source)
for(dofile in files_to_source[-dont_source]){
  source(paste0("../two_stage/vital_rate_functions/",dofile))
}
source("../buffering_experiments/buffering_fxns.R")




####
####  STATE, SPECIES, AND SIMULATION INFORMATION ----
####
# Make matrix of state-species for one-stage modeling
state_species <- data.frame(state = c(rep("idaho",4),
                                      rep("arizona",1),
                                      rep("newmexico",1)),
                            species = c(sort(c("ARTR","HECO","POSE","PSSP")),
                                        sort(c("BOER")),
                                        sort(c("SPFL"))))



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
  grow_params_adults <- readRDS(file = grow_file)[["yearonly"]]
  grow_params_kids   <- readRDS(file = grow_file)[["juveniles"]]
  
  ##  Load survival regression parameters
  surv_file <- paste0("./vital_rate_params/surv_params_",state,"_",do_species,".RDS")
  surv_params_adults <- readRDS(file = surv_file)[["yearonly"]]
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
  
  ##  Fit growth and survival models
  suppressWarnings(gdata <- fetch_grow_data(data_dir, state, do_species))
  gfit  <- lmer(logarea.t1 ~ logarea.t0 * W + (1|Group) + (1|fyear), 
                data = reformat_growth(gdata[["adults"]],state))
  suppressWarnings(sdata <- fetch_surv_data(data_dir, state, do_species))
  sfit  <- glmer(survives ~ logarea * W + (1|Group) + (1|fyear), 
                 data = reformat_survival(sdata[["adults"]],state), 
                 family = "binomial")
  
  ##  Unlist IPM specifications to named objects
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
  ####  DEFINE SIZE-YEAR EFFECTS ----
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
  
  
  ## FIGURE OUT WHERE NEW RECRUITS GO IN THE nt1 VECTOR 
  kidsize <- log(0.25)
  kidbin <- 1 + sum(u<kidsize) 
  
  
  # set to stable size distribution
  tnt <- stable_size
  # set to fix cover value
  tmp <- fix_cover*100/(h*sum(tnt*exp(u)))
  tnt <- tnt*tmp
  nt1[] <- 0
  nt1[kidbin] <- tnt[kidbin]
  nt2[] <- 0
  nt2[-kidbin] <- tnt[-kidbin]
  nt2[kidbin] <- mean(tnt[kidbin-1], tnt[kidbin+1])
  
  
  
  ####
  ####  ITERATE IPM ----
  ####
  ##  Create year vector for observed years
  the_years <- Reduce(intersect, list(rec_params$year, surv_params_adults$year_id, grow_params_adults$year_id))
  all_years <- as.numeric(the_years)
  year_save <- sort(all_years)
  tlimit <- length(year_save)+1
  
  ##  Create empty data frame with necessary dimensions to store simulation
  ipm_sims <- as.data.frame(matrix(NA, nrow = tlimit, ncol = 6)) %>%
    rename(iteration  = V1,
           year       = V2,
           cov_kids   = V3,
           cov_adults = V4,
           N_kids     = V5,
           N_adults   = V6) %>%
    mutate(iteration = 1:tlimit,
           year = c(1,year_save))
  
  ##  Plug in inital values
  ipm_sims$cov_kids[1]   <- sumCover(u,nt1,h,A) # store the % cover as cm^2/cm^2
  ipm_sims$cov_adults[1] <- sumCover(u,nt2,h,A) # store the % cover as cm^2/cm^2
  ipm_sims$N_kids[1]     <- sumN(nt1,h) # sum abundance over size stages
  ipm_sims$N_adults[1]   <- sumN(nt2,h) # sum abundance over size stages
  
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
  
  for(iiter in 2:tlimit){
    ##  Set t0 and t1 for indexing; grab the first year
    t0      <- iiter - 1
    t1      <- iiter
    do_year <- ipm_sims$year[t1] # for random year effects
    
    ##  Calculate total basal cover for recs_per_area and crowding
    cover <- fix_cover/100
    
    ##  Calculate total population size distribution
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
      Pmatrix_adults <- make.P.matrix(v=u, 
                                      Wvec=Wvec, 
                                      Gpars=grow_params_adults, 
                                      Spars = surv_params_adults, 
                                      do_year = do_year, 
                                      G=G, 
                                      S=S,
                                      surv_allsize_years = surv_size_year,
                                      grow_allsize_years = grow_size_year)
      
      ##  Project populations
      new_nt2         <- growers_nt + Pmatrix_adults%*%nt2
      new_nt1         <- rep(0,length(u)) 
      new_nt1[kidbin] <- recruits # put new recruits into their size bin, such that h*sum(new_nt1) = #recruits this year
    } # end if cover > 0
    
    ##  Update nt's for next timestep
    nt1 <- new_nt1 # kids
    nt2 <- new_nt2 # adults
    
    ##  Save population states
    ipm_sims$cov_kids[t1]   <- sumCover(u,nt1,h,A) # store the % cover as cm^2/cm^2
    ipm_sims$cov_adults[t1] <- sumCover(u,nt2,h,A) # store the % cover as cm^2/cm^2
    ipm_sims$N_kids[t1]     <- sumN(nt1,h) # sum abundance over size stages
    ipm_sims$N_adults[t1]   <- sumN(nt2,h) # sum abundance over size stages
    
    ##  Reset to low cover
    tnt <- nt1+nt2
    tmp <- fix_cover*100/(h*sum(tnt*exp(u)))
    tnt  <- tnt*tmp
    nt1[] <- 0
    nt1[kidbin] <- tnt[kidbin]
    nt2[] <- 0
    nt2[-kidbin] <- tnt[-kidbin]
    nt2[kidbin] <- mean(tnt[kidbin-1], tnt[kidbin+1])
  } # end IPM iterations
  
  
  
  ####
  ####  SAVE OUTPUT ----
  ####
  out_file <- paste0(out_dir, paste(state, do_species, isim, sep = "_"),".RDS")
  saveRDS(ipm_sims, out_file)
} # End state-species loop



