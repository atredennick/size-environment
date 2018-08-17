Author: Andrew Tredennick
Date: September 14, 2017

This directory contains R scripts for running IPMs with and without size-by-year interactions. There are two core directories, within which the file structure is similar:

1. ./two_stage/ -- contains scripts for IPM with seedling stage
2. ./one_stage/ -- contains script for IPM without seedling stage

Required packages:
install.packages(c('tidyverse','dplyr','gamlss','stats4','boot','lme4'))



OCT 6 — fit vital rate models with centered log(size). Then send SD figure to group.