################################################################################
##  fit_vital_rate_models.R -- R script that can be called from the bash script
##  "run_vital_fits.sh" to fit growth, survival, and recruitment regression
##  models for each species from our five focal data sets. To run a single
##  instance, e.g. for a single species, set 'i' manually under the "Get
##  system arguments" section below. For growth and survival, three separate
##  models are fit for each species:
##
##      1. A juvenile model for one-year-old plants/seedlings;
##      2. A (size * year) interaction  model for adult plants;
##      3. A (size + year) model for adult plants.
##
##  Models 2 and 3 allow us to compare IPM runs with and without the size by
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
args <- commandArgs(trailingOnly = F)
myargument <- args[length(args)]
myargument <- sub("-","",myargument)
i <- as.numeric(myargument)

#i <- 9 # hardset for testing

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

work_dir <- paste0(root,"/drivers/empirical/size_by_year_IPM/two_stage/runscripts/")
data_dir <- paste0(root,"/driversdata/data/")
out_dir  <- paste0(root,"/drivers/empirical/size_by_year_IPM/two_stage/vital_rate_params/")
mcmc_dir <- paste0(root,"/drivers/empirical/size_by_year_IPM/two_stage/rec_mcmc/")

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
seedling_size <- round(as.numeric(names(sort(table(exp(grow_data[["all"]]$logarea.t0)),decreasing=TRUE)[1])),5)
if(seedling_size > 0.26) seedling_size <- 0.26

grow_size_by_year <- fit_growth(df = reformat_growth(grow_data[["adults"]], state),
                                juvenile = FALSE, 
                                size_by_year = TRUE)

grow_size_and_year <- fit_growth(df = reformat_growth(grow_data[["adults"]], state),
                                 juvenile = FALSE, 
                                 size_by_year = FALSE)

### TODO: Decide on N cut off for estimate juvenile growth distribution.
if(nrow(grow_data[["juveniles"]]) > 10){
  grow_juvenile <- fit_growth(df = reformat_growth(grow_data[["juveniles"]], state),
                              juvenile = TRUE, 
                              size_by_year = FALSE, 
                              seedling_size = seedling_size)
}
if(nrow(grow_data[["juveniles"]]) < 10){
  grow_juvenile <- NULL
}



grow_sxy_params <- format_growth_params(grow_size_by_year)
grow_params     <- format_growth_params(grow_size_and_year)
growth_out_list <- list(juveniles = grow_juvenile, 
                        sizeXyear = grow_sxy_params, 
                        yearonly  = grow_params)
saveRDS(object = growth_out_list, 
        file   = paste0(out_dir,"growth_params_",state,"_",do_species,".RDS"))



####
####  FIT, FORMAT, AND SAVE SURVIVAL REGRESSIONS ----
####
surv_data <- fetch_surv_data(data_dir, state, do_species)

surv_size_by_year <- fit_survival(df = reformat_survival(surv_data[["adults"]], state),
                                  juvenile     = FALSE, 
                                  size_by_year = TRUE)

surv_size_and_year <- fit_survival(df = reformat_survival(surv_data[["adults"]], state),
                                   juvenile    = FALSE, 
                                  size_by_year = FALSE)

surv_juvenile <- fit_survival(df = reformat_survival(surv_data[["juveniles"]], state),
                              juvenile     = TRUE, 
                              size_by_year = FALSE)

surv_sxy_params      <- format_surv_params(surv_size_by_year, juveniles = FALSE, all_years = unique(surv_data$all$year))
surv_params          <- format_surv_params(surv_size_and_year, juveniles = FALSE, all_years = unique(surv_data$all$year))
surv_juvenile_params <- format_surv_params(surv_juvenile, juveniles = TRUE, all_years = unique(surv_data$all$year))
surv_out_list <- list(juveniles = surv_juvenile_params, 
                      sizeXyear = surv_sxy_params, 
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

pval_ids  <- grep("pvalue",rownames(rec_fit[[1]]))
rec_stats <- rec_fit[[1]][-pval_ids,]
recruitment_params <- format_rec_params(rec_stats, years = unique(rec_data$year))
saveRDS(object = recruitment_params,
        file   = paste0(out_dir,"rec_params_",state,"_",do_species,".RDS"))
saveRDS(object = rec_fit[[2]],
        file   = paste0(mcmc_dir,"rec_mcmc_",state,"_",do_species,".RDS"))
