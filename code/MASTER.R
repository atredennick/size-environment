# MASTER.R
#  This script is the MASTER script for reproducing the results in:
#
#                        *****************
#       Tredennick, A.T., B.J. Teller, P.B. Adler, G. Hooker, 
#     and S.P. Ellner. (2018). Size-by-environment interactions: 
#    a neglected dimension of species' responses to environmental 
#                    variation. Ecology Letters.
#                        *****************
#
#  This script sources all the other R scripts necessary to make all our
#  figures and conduct all analyses. Scripts that underly the results of 
#  particular figures are noted in each figure caption in the manuscript.
#  Note that simply running this script will take a LONG, LONG time. Thus,
#  it is better to consider this script as a recipe for the analysis, and
#  interested coders should probably just check out the sourced scripts.
#  All results files, which are generated from running these scripts, are
#  also included in this fileset. We identify the scripts that generate results
#  and those that summarize and/or plot those results.
#
# Authors:
#  Andrew Tredennick (Primary contact; atredenn@gmail.com)
#  Britta Teller (Secondary contact; brittany.j.teller@gmail.com)
#  Peter Adler
#  Giles Hooker
#  Stephen Ellner


# Load necessary packages -------------------------------------------------

needed_pckgs <- c(
  "tidyverse", "ggthemes", "ggalt", "mvtnorm", "dplyr", "lme4", "boot",
  "merTools", "RLRsim", "stringr", "cowplot", "xtable", "mgcv", "gamm4",
  "ggmcmc", "coda", "gamlss", "viridis", "stats4", "rjags"
)

# ipak function: 
#  install and load multiple R packages.
#  Check to see if packages are installed.
#  Install them if they are not, then load them into the R session.
#  https://gist.github.com/stevenworthington/3178163

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
  {
    install.packages(new.pkg,
                     dependencies = TRUE,
                     repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(sapply(pkg, require, character.only = TRUE))
}

ipak(needed_pckgs)  # run the ipak function


# 01. Create graphics for conceptual Figure 1 -----------------------------

root <- getwd()  # set root directory to this location
source("./SizeByYear_schematic.R")
rm(list = ls(all.names = TRUE))


# 02. Fit all models ------------------------------------------------------

root <- getwd()  # set root directory to this location
source("./size_by_year_models/anomaly_analysis.R")
rm(list = ls(all.names = TRUE))

# Results for Table 1 are in:
#  results/grow_model_comparisons.csv
#  results/surv_model_comparisons.csv


# 03. Plot size-by-year results (Figs. 2 & 3) -----------------------------

root <- getwd()  # set root directory to this location
source("./size_by_year_models/plot_sxy_statistics.R")
rm(list = ls(all.names = TRUE))

# Figure 2: figures/sd_anomalies_all.pdf
# Figure 3: figures/small_large_corrs.pdf

