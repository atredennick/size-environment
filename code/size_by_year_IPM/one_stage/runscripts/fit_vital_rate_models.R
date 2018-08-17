################################################################################
##  fit_vital_rate_models.R -- R script that can be called from the bash script
##  "run_vital_fits.sh" to fit growth, survival, and recruitment regression
##  models for each species from our five focal data sets. To run a single
##  instance, e.g. for a single species, set 'i' manually under the "Get
##  system arguments" section below. For growth and survival, two separate
##  models are fit for each species:
##
##      1. A (size * year) interaction  model for adult plants;
##      2. A (size + year) model for adult plants.
##
##  Models 1 and 2 allow us to compare IPM runs with and without the size by
##  year interaction. Model parameters are collated and saved in the 
##  automatically-generated folder "vital_rate_params/"; these *.RDS files are
##  later used during the IPM simulations.
##
##  ____________________________________________________________________________
##  Author: Andrew Tredennick (atredenn@gmail.com)
##  Contributors: Peter Adler, Britta Teller, Stephen Ellner
##  Date created: August 23, 2017
################################################################################


##  Clear everything...
rm(list=ls(all.names = TRUE))

####
####  GET SYSTEM ARGUMENTS ----
####
##  This is for running in batch later...set manually below for now
args <- commandArgs(trailingOnly = F)
myargument <- args[length(args)]
myargument <- sub("-","",myargument)
i <- as.numeric(myargument)

##  Set model specification manually using 'i'
##  i = 1: idaho/ARTR
##  i = 15: kansas/bohi
#i <- 1



####
####  LOAD LIBRARIES ----
####
library(gamlss)
library(stats4)
library(boot) 
library(lme4)
library(tidyverse)
library(dplyr)



####
####  SET DIRECTORIES ----
####
root <- ifelse(.Platform$OS.type=="windows","c:/repos","~/repos")
if(Sys.info()["user"] == "atredenn"){root <- "~/Repos"}

work_dir <- paste0(root,"/drivers/empirical/size_by_year_IPM/one_stage/runscripts/")
data_dir <- paste0(root,"/driversdata/data/")
out_dir  <- paste0(root,"/drivers/empirical/size_by_year_IPM/one_stage/vital_rate_params/")
mcmc_dir <- paste0(root,"/drivers/empirical/size_by_year_IPM/one_stage/rec_mcmc/")

setwd(work_dir)
dir.create(file.path(out_dir), showWarnings = FALSE)



####
####  SOURCE ALL NEEDED FUNCTIONS ----
####
files_to_source <- list.files("../vital_rate_functions/")
dont_source     <- grep("jags", files_to_source)
for(dofile in files_to_source[-dont_source]){
  source(paste0("../vital_rate_functions/",dofile))
}



####
####  STATE AND SPECIES INFORMATION ----
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

state      <- as.character(state_species[i,"state"])
do_species <- as.character(state_species[i,"species"])


####
####  FIT, FORMAT, AND SAVE GROWTH REGRESSIONS ----
####
grow_data <- fetch_grow_data(data_dir, state, do_species)

####### SPE
# e <- which(grow_data$logarea.t0 > 0); 
# grow_data_1 <- grow_data[e,]; 

grow_size_by_year <- fit_growth(df = reformat_growth(grow_data, state),
                                size_by_year = TRUE)

grow_size_and_year <- fit_growth(df = reformat_growth(grow_data, state),
                                 size_by_year = FALSE)

grow_sxy_params <- format_growth_params(grow_size_by_year)
grow_params     <- format_growth_params(grow_size_and_year)
growth_out_list <- list(sizeXyear = grow_sxy_params, 
                        yearonly  = grow_params)
saveRDS(object = growth_out_list, 
        file   = paste0(out_dir,"growth_params_",state,"_",do_species,".RDS"))


####
####  FIT, FORMAT, AND SAVE SURVIVAL REGRESSIONS ----
####

surv_data <- fetch_surv_data(data_dir, state, do_species)

###### SPE 
surv_size_by_year <- fit_survival(df = reformat_survival(surv_data, state),
                                  size_by_year = TRUE)

surv_size_and_year <- fit_survival(df = reformat_survival(surv_data, state),
                                   size_by_year = FALSE)

surv_sxy_params <- format_surv_params(surv_size_by_year, all_years = unique(surv_data$all$year))
surv_params     <- format_surv_params(surv_size_and_year, all_years = unique(surv_data$all$year))

surv_out_list <- list(sizeXyear = surv_sxy_params, 
                      yearonly  = surv_params)
saveRDS(object = surv_out_list,
        file = paste0(out_dir,"surv_params_",state,"_",do_species,".RDS"))

####
####  FIT, FORMAT, AND SAVE RECRUITMENT REGRESSIONS ----
####
rec_data <- fetch_recruit_data(data_dir, state, do_species) %>%
  arrange(year) # sort by year to keep everything in order

rec_fit <- fit_recruitment(df         = rec_data,
                           n_adapt    = 5000,
                           n_update   = 10000,
                           n_iters    = 20000,
                           n_thin     = 50,
                           model_file = "../vital_rate_functions/recruitment_jags.R")

# rec_fit <- fit_recruitment(df         = rec_data,
#                            n_adapt    = 100,
#                            n_update   = 1000,
#                            n_iters    = 1000,
#                            n_thin     = 10,
#                            model_file = "../vital_rate_functions/recruitment_jags.R")

rec_stats <- rec_fit[[1]]
recruitment_params <- format_rec_params(rec_stats, years = unique(rec_data$year))
saveRDS(object = recruitment_params,
        file   = paste0(out_dir,"rec_params_",state,"_",do_species,".RDS"))
# saveRDS(object = rec_fit[[2]],
#         file   = paste0(mcmc_dir,"rec_mcmc_",state,"_",do_species,".RDS"))
