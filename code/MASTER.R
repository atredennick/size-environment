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

source("./SizeByYear_schematic.R")


# 02. Fit all models ------------------------------------------------------


