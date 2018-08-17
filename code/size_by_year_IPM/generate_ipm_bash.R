###################################################################
##  generate_ipm_bash.R: R script to generate a bash script
##  designed to call either one-stage or two-stage IPM scripts,
##  depending on whether seedlings are different or not. The
##  bash script can then be run to simulate the correct IPM for
##  each species.
##
##  Author: Andrew Tredennick (atredenn@gmail.com)
##  Date created: 10-03-2017
###################################################################

##  Clear everything
rm(list = ls(all.names = TRUE))



####
####  LOAD LIBRARIES ----
####
library(tidyverse)
library(dplyr)
library(stringr)



####
####  SET DIRECTORIES ----
####
root     <- ifelse(.Platform$OS.type=="windows","c:/repos","~/repos")
work_dir <- paste0(root,"/drivers/empirical/size_by_year_IPM/")
setwd(work_dir) # set working directory



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



####
####  LOAD SEEDLING ANALYSIS RESULTS ----
####
suppressWarnings( # ignore message re: coerce factor to character; it's OK
  seedling_results <- readRDS("./seedling_analysis/seedling_pvalues.RDS") %>%
    dplyr::filter(vital_rate == "survival") %>%
    dplyr::mutate(ipm = ifelse(significant == "no", 
                               "./one_stage/runscripts/ipm_no_overlap_sxy.R", 
                               "./two_stage/runscripts/ipm_no_overlap_sxy.R")) %>%
    right_join(state_species, by = c("state", "species"))
)

bash_commands <- paste("Rscript", seedling_results$ipm, 1:nrow(seedling_results))
bash_script <- paste(bash_commands, collapse = "\n")
sink(file = "./run_all_ipms.sh")
cat("#!/bin/bash")
cat("\n\n")
cat(bash_script)
sink()



####
####  SENSITIVITY ANALYSIS SCRIPT ----
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

suppressWarnings( # ignore message re: coerce factor to character; it's OK
  seedling_results <- readRDS("./seedling_analysis/seedling_pvalues.RDS") %>%
    dplyr::filter(vital_rate == "survival") %>%
    dplyr::mutate(ipm = ifelse(significant == "no", 
                               "./one_stage/runscripts/sensitivity_analysis.R", 
                               "./two_stage/runscripts/sensitivity_analysis.R")) %>%
    right_join(state_species, by = c("state", "species"))
)

bash_commands <- paste("Rscript", seedling_results$ipm, 1:nrow(seedling_results))
bash_script <- paste(bash_commands, collapse = "\n")
sink(file = "./run_all_sensitivities.sh")
cat("#!/bin/bash")
cat("\n\n")
cat(bash_script)
sink()
