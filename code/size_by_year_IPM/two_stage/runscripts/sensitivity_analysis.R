################################################################################
##  ipm_no_overlap_sxy.R: R script that can be called from bash script
##  "run_ipms.sh" to run single species Integral Projection Model (IPM) 
##  simulations for each of our 15 focal species. IPMs are run with and without 
##  size-by-year interactions to analyze the effect of size-specific climate 
##  responses on long-term population dynamics. The code returns and saves a 
##  data frame of species cover over time, for each species and each model type 
##  (interaction or not). This code requires several utility functions from the 
##  './ipm_functions/' folder in this directory.
##
##  ____________________________________________________________________________
##  Author: Andrew Tredennick (atredenn@gmail.com)
##  Contributors: Peter Adler, Britta Teller, Stephen Ellner
##  Date created: August 28, 2017
################################################################################


##  Clear everything
rm(list = ls(all.names = TRUE))

####
####  GET SYSTEM ARGUMENTS ----
####
args       <- commandArgs(trailingOnly = F)
myargument <- args[length(args)]
myargument <- sub("-","",myargument)
i          <- as.numeric(myargument)
# i <- 4 # hard set for testing



####
####  LOAD LIBRARIES ----
####
library(tidyverse)
library(dplyr)
library(gamlss)
library(stats4)
library(boot) 
library(lme4)
library(directlabels)



####
####  SET DIRECTORIES ----
####
root     <- ifelse(.Platform$OS.type=="windows","c:/repos","~/repos")
work_dir <- paste0(root,"/drivers/empirical/size_by_year_IPM/two_stage/")
data_dir <- paste0(root,"/driversdata/data/")
out_dir  <- paste0(root,"/drivers/empirical/size_by_year_IPM/ipm_results/sensitivities/")

setwd(work_dir) # set working directory
dir.create(file.path(out_dir), showWarnings = FALSE) # create a folder for results



####
####  SOURCE ALL NEEDED FUNCTIONS ----
####
files_to_source <- list.files("./ipm_functions/")
for(dofile in files_to_source){
  source(paste0("./ipm_functions/",dofile))
}



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

## Read in size-by-year statistical results to pull right model
grow_stats <- read.csv("../../size_by_year_models/results/grow_model_comparisons.csv") %>%
  mutate(
    growmod = ifelse(obs_lrt > uquant, "sizeXyear", "yearonly"),
    state = tolower(state),
    species = as.character(species)
  ) %>%
  dplyr::select(state, species, growmod)

surv_stats <- read.csv("../../size_by_year_models/results/surv_model_comparisons.csv") %>% 
  mutate(
    survmod = ifelse(obs_lrt > uquant, "sizeXyear", "yearonly"),
    state = tolower(state),
    species = as.character(species)
  ) %>%
  dplyr::select(state, species, survmod)

state_species <- state_species %>% 
  mutate(
    state = as.character(state),
    species = as.character(species)
  ) %>% 
  left_join(grow_stats, by = c("state", "species")) %>% 
  left_join(surv_stats, by = c("state", "species"))

state        <- as.character(state_species[i,"state"])
do_species   <- as.character(state_species[i,"species"])
grow_size_by_year <- state_species[i,"growmod"]
surv_size_by_year <- state_species[i,"survmod"]



####
####  LOAD VITAL RATE PARAMETERS ----
####
##  Load growth regression parameters
grow_file <- paste0("./vital_rate_params/growth_params_",state,"_",do_species,".RDS")
grow_params_adults <- readRDS(file = grow_file)[[grow_size_by_year]]
grow_params_kids   <- readRDS(file = grow_file)[["juveniles"]]

##  Load survival regression parameters
surv_file <- paste0("./vital_rate_params/surv_params_",state,"_",do_species,".RDS")
surv_params_adults <- readRDS(file = surv_file)[[surv_size_by_year]]
surv_params_kids   <- readRDS(file = surv_file)[["juveniles"]]

##  Load recruitment regression parameters
rec_file <- paste0("./vital_rate_params/rec_params_",state,"_",do_species,".RDS")
rec_params <- readRDS(file = rec_file)
rec_size_params <- get_rec_size(data_dir,state,do_species)


####
####  SET UP IPM SPECIFICATIONS AND INITIAL CONDITIONS ----
####
if(state %in% c("arizona", "kansas", "newmexico")){
  grow_data <- read.csv(paste0(data_dir,state,"/speciesData/",do_species,"/growDnoNA_reduced.csv"))
  surv_data <- read.csv(paste0(data_dir,state,"/speciesData/",do_species,"/survD_reduced.csv"))
}
if(state %in% c("idaho", "montana")){
  grow_data <- read.csv(paste0(data_dir,state,"/speciesData/",do_species,"/growDnoNA.csv"))
  surv_data <- read.csv(paste0(data_dir,state,"/speciesData/",do_species,"/survD.csv"))
}

spec_list  <- get_ipm_specs(grow_data)
inits_list <- get_initial_conditions(tlimit = spec_list[["tlimit"]],
                                     u = spec_list[["u"]],
                                     h = spec_list[["h"]],
                                     A = spec_list[["A"]],
                                     surv_data = surv_data)

##  Unlist IPM specifications and initial conditions to named objects
for(i in 1:length(names(spec_list))){
  assign(names(spec_list)[[i]], spec_list[[i]])
}
for(i in 1:length(names(inits_list))){
  assign(names(inits_list)[[i]], inits_list[[i]])
}
 
##  Hard set -- remember to remove!!!
# tlimit <- 400
# burn_in <- 50



####
####  SETUP IPM ----
####
##  Sample years
the_years <- Reduce(intersect, list(as.numeric(rec_params$year), 
                                    as.numeric(surv_params_adults$year_id), 
                                    as.numeric(grow_params_adults$year_id)))
all_years <- as.numeric(the_years)
year_save <- sample(all_years, tlimit+1, replace = TRUE)

##  Plug in initial values
cov_kids <- cov_save_kids[1]
cov_adults <- cov_save_adults[1]

##  Calculate crowding variables
Wvars <- get_crowd_vars(data_dir = data_dir, 
                        state    = state, 
                        species  = do_species, 
                        b.r      = b.r)

##  Calculate the seedling growth size distribution; this is constant
kidsize <- log(0.25)
kidbin <- 1 + sum(u<kidsize) 
grow_dist  <- get_growth_kids(u,u,grow_params_kids,kidsize)
if(round(sum(h*grow_dist))!=1){
  stop("Check grow_dist function, h*grow_dist should equal 1","\n")
}



####
####  SENSITIVITY ANALYSIS STORAGE VECTORS AND MATRICES ----
####
##  Chop off one year from tlimit to set size for simulations
n_est <- tlimit

##  A vector to hold the growth rates calculated at each time step
rt.N <- rep(NA,n_est)

##  Specify the initial population structure vector
wvec <- rep(1,bigM)
wvec <- wvec/bigM

##  Specify the initial reproductive value vector
vvec <- rep(1, bigM)
vvec <- vvec/bigM



####
####  ITERATE IPM ----
####
##  Get wt time series and growth rate
pb <- txtProgressBar(min=2, max=tlimit, char="+", style=3, width=65)

### Get wt time series and growth rate ###
wt <- matrix(1/bigM, nrow = n_est+1, ncol = bigM+1) # add one to bigM for seedling class 
all_Ks <- list()

for(iiter in 1:(n_est+1)){
  do_year <- year_save[iiter] # for random year effects
  
  ##  Calculate total population size distribution
  tnt  <- nt1+nt2
  expu <- exp(u)
  
  cover <- sum(cov_kids, cov_adults)
  
  ##  Calculate crowding vector across sizes
  if(do_species %in% c("BOGR","SCSC")){
    Cr <- approxfun(b.r,h*c(0,cumsum(expu*tnt)),rule=2) # avoids integration issues
  }else{
    Cr <- splinefun(b.r,h*c(0,cumsum(expu*tnt)),method="natural")
  }
  
  Wvec <- makeW(Wvars,u,A,Cr)

  if(cover > 0){ # Don't do anything if species is extinct, for speed
    ##  Calculate recruits
    rpa      <- get_rec_number(params = rec_params, cover = cover, do_year = do_year)
    recruits <- f(v=u,u=u,Rpars=rec_size_params,rpa = rpa) 
    
 
    ##  Calculate growth*survival kernel / projection matrix
    # Kids/seedlings
    
    ############# Possible problem ######################################################################
    ## SPE is worried about the next line. get_survival_kids() should be one number, because seedlings are
    ## all the same size. But get_survival_kids is returning a vector. I think that what should be happening
    ## here is that the right element is pulled out of get_survival_kids, for the survival resulting
    ## from the W value for by the size bin that kids fall into. 
    #####################################################################################################
    growers_nt <- grow_dist * sumN(nt1,h) * get_survival_kids(u,u,Wvec,surv_params_adults,do_year)
    
    ### In other words, I think it should be 
    kid_survival_vec = get_survival_kids(u,u,Wvec,surv_params_adults,do_year)
    growers_nt <- grow_dist * sumN(nt1,h) * kid_survival_vec[kidbin]
    ### END SPE 
    
        # Adult plants
    Pmatrix_adults <- make.P.matrix(v=u, Wvec=Wvec, Gpars=grow_params_adults, Spars = surv_params_adults, do_year = do_year, G=G, S=S)
    
    ##  Project populations
    new_nt2         <- growers_nt + Pmatrix_adults%*%nt2
    new_nt1         <- rep(0,length(u)) 
    new_nt1[kidbin] <- recruits # put new recruits into their size bin, such that h*sum(new_nt1) = #recruits this year
  } # end if cover > 0
  
  if(iiter != (n_est+1)){
    ##  Make expanded K

    ##### SPE: calculate per-capita fecundity (in recruits), by size class of adults 
    ## These are the numbers that, multiplied by the nt2 vector, sum to the total number of new recruits 
    ## Note: this makes the 'recruits' category the NUMBER of recruits, not the density per size category. 
    fec_by_class = recruits*exp(u)/sum(exp(u)*nt2)
    if(iiter%%100==1) cat(iiter,"  ",sum(fec_by_class*nt2)-recruits," should equal 0","\n") 
  
    ## SPE: Now put capita fecundities into the top of K
    seedling_row <- matrix(fec_by_class, nrow = 1, ncol = ncol(Pmatrix_adults))
    
    # Create empty column of zeros for left of K (all seedlings)  
    seedling_col <- matrix(0, nrow = (ncol(Pmatrix_adults) + 1), ncol = 1)
        
    ### SPE: this year's seedlings move up to be adults.  
    ### sum(h*grow_dist)=1, so this distributes surviving seedlings as a size distribution, preserving 
    ### their total number when the distribution is integrated by midpoint rule. 
    ### The first entry of seedling_col is zero - seedlings don't remain seedlings.
    ###
    ### Note that this includes my proposed correction(?) to the survival of seedlings. 
    seedling_col[2:nrow(seedling_col),1] <- grow_dist*kid_survival_vec[kidbin]
    
    K <- rbind(seedling_row, Pmatrix_adults) # add the new row to top of Pmatrix
    K <- cbind(seedling_col, K) # add the new column to the left end of Pmatrix
    ######### END SPE 
 
    #matrix.image(K)
    
    ## calculating the population size vector, one-step growth rates, and the stage distribution vectors.
    nprime <- new_nt1 + new_nt2 # calculate the next population size vector.
    rt.N[iiter] <- sum(nprime)/sum(tnt) # calculate the one-step growth rate.
    wt[iiter+1,] <- K%*%wt[iiter,]  # calculate the next population structure vector.
    wt[iiter+1,] <- wt[iiter+1,] / sum(wt[iiter+1,]) # scale the population structure vector.
    all_Ks[[iiter]] <- K # save the kernel for year iiter
    
    ##  Update nt's for next timestep
    nt1 <- new_nt1 # kids
    nt2 <- new_nt2 # adults
    
    ##  Save population states
    cov_kids   <- sumCover(u,nt1,h,A) # store the % cover as cm^2/cm^2
    cov_adults <- sumCover(u,nt2,h,A) # store the % cover as cm^2/cm^2
  }
  
  if(iiter == (n_est+1)){
    all_Ks[[iiter]] <- K # save the kernel for year iiter
  }
  
  setTxtProgressBar(pb, iiter)
} # end IPM iterations


### Get vt time series ###
vt <- matrix(1/bigM, nrow = n_est+1, ncol = bigM+1) 
for(iiter in (n_est+1):2){
  K <- all_Ks[[iiter]]
  vt[iiter-1,] <- vt[iiter,] %*% K
  vt[iiter-1,] <- vt[iiter-1,] / sum(vt[iiter-1,])
  if (iiter%%100 == 0) cat("vt",iiter,"\n")
} # end IPM iterations


### Calculate stochastic sensitivities ###
sens.s <- matrix(0,bigM+1,bigM+1)
elas.s <- matrix(0,bigM+1,bigM+1)
for(iiter in burn_in:(n_est - burn_in)){
  K <- all_Ks[[iiter]]
  vt1.wt <- vt[iiter+1,] %*% t(wt[iiter,])
  vt1.K.wt <- sum(vt[iiter+1,] * (K %*% wt[iiter,]))
  sens.s <- sens.s + vt1.wt / vt1.K.wt
  elas.s <- elas.s + K * (vt1.wt / vt1.K.wt)
}

loglambda_s <- mean(log(rt.N[burn_in:tlimit]))
lambda_s <- exp(loglambda_s)
sens.s <- lambda_s * sens.s / (n_est - 2 * burn_in + 1)
elas.s <- elas.s / (n_est - 2 * burn_in + 1)

output <- list(sensitivity = sens.s, elasticity = elas.s, size_bins = u, wt = wt, vt = vt)
saveRDS(output, file = paste0(out_dir,state,"_",do_species,"_sens.RDS"))

### Plot the sensitivity surface ###
# matrix.image <- function(A, x=NULL, y=NULL, col=rainbow(100,start=0.67,end=0),
#                          bw=FALSE, do.contour=FALSE, do.legend=TRUE,...) {
#   if(do.legend) layout(mat=cbind(matrix(1,5,5),rep(2,5)));
#   par(mar=c(6,5,3,2));
#   if(is.null(x)) x=1:ncol(A);
#   if(is.null(y)) y=1:nrow(A);
#   nx=length(x); ny=length(y);
#   x1=c(1.5*x[1]-0.5*x[2],1.5*x[nx]-0.5*x[nx-1]);
#   y1=c(1.5*y[1]-0.5*y[2],1.5*y[ny]-0.5*y[ny-1]);
#   if(bw) col=grey( (200:50)/200 );
#   image(list(x=x,y=y,z=A),xlim=x1,ylim=rev(y1),col=col,cex.axis=1,cex.lab=1,bty="u",...);
#   abline(v=range(x1)); abline(h=range(y1));
#   if(do.contour) {
#     labs <- round(tapply(as.numeric(A), cut(as.numeric(A), 5), max),1)
#     contour(x,y,A,nlevels=5,labcex=0.5,add=TRUE);
#   }
#   
#   if(do.legend) {
#     l.y=seq(min(A),max(A),length=100);
#     par(mar=c(6,2,3,1))
#     image(list(x=1:2,y=l.y,z=rbind(l.y,l.y)),col=col,bty="o",xaxt="n",yaxt="n");
#     axis(side=2,cex.axis=1,at=pretty(seq(min(A),max(A),length=10)));
#   }
# }
# par(mfrow=c(1,2),tcl=-0.2,mgp=c(2,0.5,0),mar=c(3,4,1,1),cex.lab=1.2,oma=c(0, 0, 0, 0), las=1)
# matrix.image(sqrt(t(sens.s[2:101,2:101])),do.contour = TRUE, bw=TRUE, 
#              xlab = "Size (t), z", ylab = "Size (t+1), z\'", do.legend = FALSE)
# matrix.image(t(elas.s[2:101,2:101]),do.contour = TRUE, bw=TRUE, 
#              xlab = "Size (t), z", ylab = "Size (t+1), z\'", do.legend = FALSE)
