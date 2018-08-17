################################################################################
##  rec_mcmc_diagnostics.R: script to check MCMC chains for stationarity
##  and convergence
##
##  Author: Andrew Tredennick (atredenn@gmail.com)
##  Date created: 10-31-2017
################################################################################

##  Clear everything
rm(list = ls(all.names = TRUE))


####
####  LOAD LIBRARIES ----
####
library(tidyverse) # Data manipulation
library(stringr)   # Workging with strings
library(ggthemes)  # Pleasing ggplot2 themes
library(ggmcmc)     # Some more geoms for ggplot2
library(cowplot)   # Combining ggplot2 objects
library(coda)      # MCMC diagnostics



####
####  SET DIRECTORIES ----
####
root        <- ifelse(.Platform$OS.type=="windows","c:/repos","~/repos")
work_dir    <- paste0(root,"/drivers/empirical/size_by_year_IPM/")
loc_fig_dir <- paste0(root,"/drivers/empirical/size_by_year_IPM/analysis/figures/")


setwd(work_dir) # set working directory
dir.create(file.path(loc_fig_dir), showWarnings = FALSE) # create a folder for figures




####
####  READ IN MCMC RESULTS AND CHECK ----
####
file_paths <- paste0("./one_stage/rec_mcmc/",list.files("./one_stage/rec_mcmc/"))
file_paths <- c(file_paths, paste0("./two_stage/rec_mcmc/",list.files("./two_stage/rec_mcmc/")))
if(length(file_paths) != 15){
  stop("Missing some species' MCMC chains.")
}

gelman_out <- {}
for(do_file in file_paths){
  ##  File information
  state_spp <- str_sub(do_file, 31, end = nchar(do_file))
  state_spp <- str_sub(state_spp, -nchar(state_spp),-5)
  state <- str_split(state_spp,"_")[[1]][1]
  species <- str_split(state_spp,"_")[[1]][2]
  
  ##  Caclulate multivariate Rhat and save traceplots
  tmp_mcmc        <- readRDS(do_file)
  gelman_multivar <- gelman.diag(tmp_mcmc)[[2]]
  pdf(paste0(loc_fig_dir, state,"_",species,"_mcmc.pdf"))
  plot(tmp_mcmc, density = FALSE)
  dev.off()
  
  ##  Calculate P-values for posterior predictive checks
  mcmc_stats <- summary(tmp_mcmc)$stat
  
  ##  Store output
  tmp_out <- data.frame(state = state,
                        species = species,
                        gelman = gelman_multivar,
                        pvalue_mean = mcmc_stats[which(rownames(mcmc_stats)=="pvalue_mean"),"Mean"],
                        pvalue_sd = mcmc_stats[which(rownames(mcmc_stats)=="pvalue_sd"),"Mean"])
  gelman_out <- rbind(gelman_out, tmp_out)
}



