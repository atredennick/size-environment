Author: Andrew Tredennick
Date: September 14, 2017

This directory contains R scripts for fitting vital rate regressions and running an IPM with no seedling stage. Vital rate regressions are fit with and without size-by-year interactions. The IPM is run similarly, with and without the size-by-year interactions. The main scripts are in the ./runscripts/ directory. File structure is described below.

1. ./ipm_functions/ — Contains functions for estimating crowding in the IPM; setting up IPM parameters (e.g., nt vector); predicting growth, survival, and recruitment; and calculating the IPM iteration matrix (Kmatrix).

2. ./ipm_results/ — This directory is made on the fly to store IPM simulation results.

3. ./runscript/ — Contain scripts for (i) fitting all vital rate regressions with and without size-by-year effects and (ii) running the IPM. The file names should make this obvious. There is also a bach script (run_vital_fits.sh) that calls runs a batch of vital rate fitting R scripts for each species (avoids looping in the R script).

4. ./vital_rate_functions/ — This directory contains scripts with functions for importing and formatting demographic data, fitting statistical models, and formatting model fits for import to IPMs.

5. ./vital_rate_params/ — This directory is made on the fly to store vital rate parameters as *.RDS files.


!!!!!!!
NOTES:
!!!!!!!
1. Divergent integral errors for Montana::BOGR and Kansas::BOHI; related to making W.


