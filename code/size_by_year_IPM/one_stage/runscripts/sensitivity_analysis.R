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
args <- commandArgs(trailingOnly = F)
myargument <- args[length(args)]
myargument <- sub("-","",myargument)
i <- as.numeric(myargument)
# i <- 21 # hard set for testing; Montana::HECO


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
out_dir  <- paste0(root,"/drivers/empirical/size_by_year_IPM/ipm_results/sensitivities/")

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

## Read in size-by-year statistical results to pull right model
grow_stats <- read.csv("../../size_by_year_models/results/grow_model_comparisons.csv") %>%
  mutate(
    growmod = ifelse(obs_lrt > uquant, "sizeXyear", "yearonly"),
    state = tolower(state),
    species = as.character(species)
  ) %>%
  dplyr::select(state, species, growmod)

surv_stats <- read.csv("../../size_by_year_models/results/surv_model_comparisons.csv") %>% 
  mutate(
    survmod = ifelse(obs_lrt > uquant, "sizeXyear", "yearonly"),
    state = tolower(state),
    species = as.character(species)
  ) %>%
  dplyr::select(state, species, survmod)

state_species <- state_species %>% 
  mutate(
    state = as.character(state),
    species = as.character(species)
  ) %>% 
  left_join(grow_stats, by = c("state", "species")) %>% 
  left_join(surv_stats, by = c("state", "species"))

state        <- as.character(state_species[i,"state"])
do_species   <- as.character(state_species[i,"species"])
grow_size_by_year <- state_species[i,"growmod"]
surv_size_by_year <- state_species[i,"survmod"]



####
####  LOAD VITAL RATE PARAMETERS ----
####
##  Load growth regression parameters
grow_file <- paste0("./vital_rate_params/growth_params_",state,"_",do_species,".RDS")
grow_params <- readRDS(file = grow_file)[[grow_size_by_year]]

##  Load survival regression parameters
surv_file <- paste0("./vital_rate_params/surv_params_",state,"_",do_species,".RDS")
surv_params <- readRDS(file = surv_file)[[surv_size_by_year]]

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

##  Hard set -- remember to remove!!!
# tlimit <- 500
# burn_in <- 50

####
####  SET UP IPM ----
####
##  Sample years
the_years <- Reduce(intersect, list(rec_params$year, surv_params$year_id, grow_params$year_id))
all_years <- as.numeric(the_years)
year_save <- sample(all_years, tlimit+1, replace = TRUE)

##  Plug in initial values
cover <- sumCover(u,nt,h,A)

##  Calculate crowding variables
Wvars <- get_crowd_vars(data_dir = data_dir, 
                        state    = state, 
                        species  = do_species, 
                        b.r      = b.r)

max_size <- log(maxSize)
fecundity_cutoff <- sum(u<max_size) 

####
####  SENSITIVITY ANALYSIS STORAGE VECTORS AND MATRICES ----
####
##  Chop off one year from tlimit to set size for simulations,
##  becuase simulations start in year 2
n_est <- tlimit

##  A vector to hold the growth rates calculated at each time step
rt.N <- rep(NA,n_est)

##  Specify the initial population structure vector
wvec <- rep(1,bigM)
wvec <- wvec/bigM

##  Specify the initial reproductive value vector
vvec <- rep(1, bigM)
vvec <- vvec/bigM

####
####  ITERATE IPM ----
####
##  Run IPM from t = 2 to t = tlimit, referenced by 'iiter'

meanKmatrix <- matrix(0,bigM,bigM); 

pb <- txtProgressBar(min=2, max=tlimit, char="+", style=3, width=65) # progress bar

### Get wt time series and growth rate ###
wt <- matrix(1/bigM, nrow = n_est+1, ncol = bigM) 
all_Ks <- list()
for(iiter in 1:(n_est+1)){
  do_year <- year_save[iiter] # for random year effects

  ##  Calculate total population size
  N    <- sumN(nt,h)
  expu <- exp(u)
  
  ##  Calculate crowding vector across sizes
  if(do_species %in% c("BOGR","SCSC")){
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
                             fecundity_cutoff = fecundity_cutoff)	
    
    meanKmatrix = meanKmatrix+Kmatrix; 
    
    ##  Project population
    nprime <- Kmatrix%*%nt
    
  } # end if cover > 0
  
  if(iiter != (n_est+1)){
    ## calculating one-step growth rates, and the wt vector
    rt.N[iiter] <- sum(nprime)/sum(nt) # calculate the one-step growth rate.
    wt[iiter+1,] <- Kmatrix%*%wt[iiter,]  # calculate the next population structure vector.
    wt[iiter+1,] <- wt[iiter+1,] / sum(wt[iiter+1,])
    all_Ks[[iiter]] <- Kmatrix # save the kernel for year iiter
    
    ##  Update nt's for next timestep
    nt <- nprime
    
    ##  Update cover
    cover <- sumCover(u,nt,h,A)
  }
  
  if(iiter == (n_est+1)){
    all_Ks[[iiter]] <- Kmatrix # save the kernel for year iiter
  }
  
  ##  Update progess
  setTxtProgressBar(pb, iiter)
} # end IPM iterations
meanKmatrix = meanKmatrix/iiter; 

### Get vt time series ###
vt <- matrix(1/bigM, nrow = n_est+1, ncol = bigM) 
for(iiter in (n_est+1):2){
  K <- all_Ks[[iiter]]
  vt[iiter-1,] <- vt[iiter,] %*% K
  vt[iiter-1,] <- vt[iiter-1,] / sum(vt[iiter-1,])
  if (iiter%%100 == 0) cat("vt",iiter,"\n")
} # end IPM iterations


### Calculate stochastic sensitivities ###
sens.s <- matrix(0,bigM,bigM)
elas.s <- matrix(0,bigM,bigM)
for(iiter in burn_in:(n_est - burn_in)){
  K <- all_Ks[[iiter]]
  vt1.wt <- vt[iiter+1,] %*% t(wt[iiter,])
  vt1.K.wt <- sum(vt[iiter+1,] * (K %*% wt[iiter,]))
  sens.s <- sens.s + vt1.wt / vt1.K.wt
  elas.s <- elas.s + K * (vt1.wt / vt1.K.wt)
}

loglambda_s <- mean(log(rt.N[burn_in:tlimit]))
lambda_s <- exp(loglambda_s)
sens.s <- lambda_s * sens.s / (n_est - 2 * burn_in + 1)
elas.s <- elas.s / (n_est - 2 * burn_in + 1)

####
####  SAVE THE OUTPUT FOR PLOTTING ----
####
output <- list(sensitivity = sens.s, elasticity = elas.s, size_bins = u, wt = wt, vt = vt)
saveRDS(output, file = paste0(out_dir,state,"_",do_species,"_sens.RDS"))
# output<- readRDS(file = paste0(out_dir,state,"_",do_species,"_sens.RDS"))


####
####  SAVE THE OUTPUT FOR PLOTTING ----
####

####### Image plot of sqrt(elasticity) with contours on arithmetic scale. 
####### Note: the loaded matrix.image has the t() built in, so don't also do it here.

#graphics.off(); 
# matrix.image(sqrt(elas.s), u, u, do.contour = FALSE, bw=TRUE,
#               xlab = "Size (t), z", ylab = "Size (t+1), z'", do.legend = FALSE)
# contour(u,u,t(elas.s),nlevels=8,labcex=1.2,add=TRUE);                 
# title(paste(state,do_species,sep="::"))
# 
# # look at the mean K matrix to make sure the plotting is rightside-up 
# #dev.new();
# matrix.image(log10(meanKmatrix), u, u, do.contour = TRUE, bw=TRUE,
#               xlab = "Size (t), z", ylab = "Size (t+1), z'", do.legend = FALSE)
# abline(0,1,col="red",lty=2,lwd=2); 
# title(paste(state,do_species,"log10Mean Kernel",sep="::"))
