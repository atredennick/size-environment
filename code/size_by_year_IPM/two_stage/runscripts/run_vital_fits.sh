#!/bin/bash 

## Idaho species
Rscript fit_vital_rate_models.R 1
Rscript fit_vital_rate_models.R 2
Rscript fit_vital_rate_models.R 3
Rscript fit_vital_rate_models.R 4

## Arizona BOER
Rscript fit_vital_rate_models.R 9

## New Mexico SPFL
Rscript fit_vital_rate_models.R 12
