####################################################################
##  collate_seedling_fits.R: R script to extract and summarize
##  the statistical fits for survival and growth when we include
##  "seedling" as a covariate. A significant "seedling" effect
##  indicates seedlings are different we should mode them as such.
##
##  Author: Andrew Tredennick (atredenn@gmail.com)
##  Date created: 9-28-2017
####################################################################


####
####  PRELIMINARIES ----
####

## Provide directory information
setwd(paste0(root,"/size_by_year_IPM/seedling_analysis/")); 

##  Load packages
library(tidyverse)
library(dplyr)
library(stringr)
library(mgcv)
library(lme4)



####
####  LOAD THE FITTED MODELS ----
####
if (file.exists("seedling_pvalues.RDS")) file.remove("seedling_pvalues.RDS")

rds_files <- as.data.frame(list.files()[grep("*.RDS",list.files())]) %>%
  dplyr::rename("fname" = !!names(.[1])) %>%
  separate(fname, c("tmp1","tmp2","tmp3","state"), sep="_") %>%
  separate(state, c("state_name","ext"), sep="[.]")

states <- rds_files$state_name
rds_files <- list.files()[grep("*.RDS",list.files())]
out_df <- {} # empty object for storage
for(do_state in states){
  state_fits <- readRDS(rds_files[grep(do_state,rds_files)])
  surv_fits  <- state_fits[["Sfits"]]
  grow_fits  <- state_fits[["Gfits"]]
  species    <- names(surv_fits)
  for(do_species in species){
    species_grow <- grow_fits[[do_species]]
    grp_rms <- grep("Group",rownames(species_grow))
    species_grow <- species_grow[-grp_rms,]
    grow_pval <- species_grow[grep("seedling",rownames(species_grow)),"Pr(>|t|)"]
    
    species_surv <- surv_fits[[do_species]]
    grp_rms <- grep("Group",rownames(species_surv))
    species_surv <- species_surv[-grp_rms,]
    surv_pval <- species_surv[grep("seedling",rownames(species_surv)),"Pr(>|z|)"]
    
    if(length(grow_pval)<1){ grow_pval <- NA }
    if(length(surv_pval)<1){ surv_pval <- NA }
    tmpdf <- data.frame(state = do_state,
                        species = do_species,
                        vital_rate = c("growth","survival"),
                        seedling_pval = c(grow_pval,surv_pval))
    out_df <- rbind(out_df,tmpdf)
  }
}

sig_ones <- which(out_df$seedling_pval < 0.05)
out_df$significant <- "no"
out_df$significant[sig_ones] <- "yes"

saveRDS(out_df,"./seedling_pvalues.RDS")

setwd("../../")

