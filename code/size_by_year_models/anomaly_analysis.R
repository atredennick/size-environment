################################################################################
##  anomaly_analysis.R: R script to fit random intercept and random intercept/
##  slope models for growth and survival, then uses the models to calculate
##  yearly anomalies for large and small plants. These anomalies define the
##  impact of the random year effects, and help distinguish where size-dependent
##  environmental responses are operating.
##
## -----------------------------------------------------------------------------
##  Created by: Andrew Tredennick (atredenn@gmail.com)
##              Brittany Teller
##              Stephen Ellner
##              Giles Hooker
##
##  Date created: December 5, 2017
################################################################################



####
####  LOAD LIBRARIES -----------------------------------------------------------
####
library(tidyverse) # Data science packages
library(dplyr)     # Data summarizing and manipulation functions
library(lme4)      # Fitting mixed-effects models
library(boot)      # Bootstrapping capabilities
library(merTools)  # Post-hoc analyses of mixed-effects models
library(RLRsim)    # Testing mixed-effects models



####
####  SET DIRECTORIES ----------------------------------------------------------
####
mainDir   <- root
utilDir   <- paste0(mainDir, "/size_by_year_models/utilities")
datadrivs <- "../data/"
outDir    <- paste0(mainDir, "../results/")

source(paste(utilDir,"/fetchDemoData.R", sep="")) # fetch data function
source(paste0(utilDir,"/TrimQuads_BT.R")) # Distinguish seedlings (rising 2yo) function


####
####  DEFINE STATE AND SPECIES DATAFRAME ---------------------------------------
####
states  <- c(rep("Idaho",4),
             rep("Montana",4),
             rep("Arizona",2),
             rep("NewMexico",2),
             rep("Kansas",3))
species <- c("ARTR","HECO","POSE","PSSP",
             "BOGR","HECO","PASM","POSE",
             "BOER","BORO",
             "BOER","SPFL",
             "ANGE","BOCU","BOHI")
state_species <- data.frame(state = states, species = species)



####
####  FIT GROWTH MODELS AND CALCULATE ANOMALIES --------------------------------
####
##  Define empty storage objects
rye_year_anoms <- {}
rye_sxy_anoms  <- {}
mod_comp_all   <- {}

## Loop over state-species
for(i in 1:nrow(state_species)){
  print(paste("Working on", as.character(state_species[i, 1]), as.character(state_species[i, 2])))
  ##  Set data directory
  do_state   <- state_species[i, "state"]
  do_species <- state_species[i, "species"]
  data_dir   <- paste0(datadrivs, do_state)
  
  ## Import demographic data
  dat <- fetchDemoData(do_species,
                       do_state,
                       grow=TRUE,
                       surv=FALSE,
                       recr=FALSE,
                       data_dir)
  gdata <- dat$data %>%
    dplyr::select(ID, year, logarea.t1, logarea.t0, W, age, Group, quad) %>%
    mutate(year = as.factor(year))
  if(min(gdata$W) == 0){gdata$W <- gdata$W+.0001}
  gdata <- gdata[gdata$logarea.t1 > log(0.25), ] # get rid of all possible seedlings
  
  ##  Drop seedlings and likely seedlings
  gdata <- trimQuadrats(data=gdata, 
                        skip.years = TRUE, 
                        df="df23", 
                        trans="grow")$data
  
  ##  Fit models
  mod_rye_year <- lmer(logarea.t1~logarea.t0*W+(1|Group)+(1|year), data = gdata)
  mod_rye_size_year <- lmer(logarea.t1~logarea.t0*W+(1|Group)+(logarea.t0|year), data = gdata)
  
  ## Loop over many iterations, simulate, compute LL: generates null distribution
  num_iters <- 500
  lrt <- numeric(num_iters)
  
  pb <- txtProgressBar(min = 1, max = num_iters, style = 3)
  for(iiter in 1:num_iters){
    gtmp <- gdata # make copy of gdata
    y_star <- simulate(mod_rye_year) # simulate model
    gtmp$logarea.t1 <- as.numeric(y_star[,1])
    
    # Fit models with simulated data from the null hypothesis
    tmp_null <- lmer(logarea.t1~logarea.t0*W+(1|Group)+(1|year), data = gtmp)
    tmp_full <- lmer(logarea.t1~logarea.t0*W+(1|Group)+(logarea.t0|year), data = gtmp)
    
    # Calculate log likelihoods
    ll_null <- logLik(tmp_null)
    ll_full <- logLik(tmp_full)
    lrt[iiter] <- 2*as.numeric(ll_full - ll_null)
    setTxtProgressBar(pb, iiter)
  }
  
  ## Calculate upper 95th quantile as value to exceed, for one-sided test
  upper <- as.numeric(quantile(lrt, probs = 0.95))
  
  ## Calculuate observed log likelihoods
  ll_null <- logLik(mod_rye_year)
  ll_full <- logLik(mod_rye_size_year)
  obs_lrt <- 2*as.numeric(ll_full - ll_null)
  
  ##  Store likelihood ratio tests
  mod_comp_tmp <- data.frame(state = do_state,
                             species = do_species,
                             obs_lrt = obs_lrt,
                             uquant = upper)
  mod_comp_all <- rbind(mod_comp_all, mod_comp_tmp)
  
  
  
  ####
  ####  Anomalies for random year on intercept model
  ####
  ##  Bootstrap effects
  fixed   <- FEsim(mod_rye_year, n.sims = 500) # bootstrapped fixed effects
  randoms <- REsim(mod_rye_year, n.sims = 500) # bootstrapped random effects

  ##  Calculate 10th and 90th percentiles of plant size
  small <- as.numeric(quantile(gdata$logarea.t0, probs=c(0.1)))
  large <- as.numeric(quantile(gdata$logarea.t0, probs=c(0.9)))

  ##  Calculate average crowding for large and small plants
  Wvals    <- c(NA,NA)
  Wvals[1] <- mean(gdata$W[which(gdata$logarea.t0 < small)])
  Wvals[2] <- mean(gdata$W[which(gdata$logarea.t0 > large)])

  ##  Calculate main effects for large and small plants
  mainSmall <- fixed$mean[1] + small*fixed$mean[2] + Wvals[1]*fixed$mean[3] + small*Wvals[1]*fixed$mean[4]
  mainLarge <- fixed$mean[1] + large*fixed$mean[2] + Wvals[2]*fixed$mean[3] + large*Wvals[2]*fixed$mean[4]

  ##  Extract random intercept parameters
  int <- randoms[randoms$term=="(Intercept)" & randoms$groupFctr=="year",]
  int$meanint <- int$mean
  int$sdint   <- int$sd

  ##  Calculate the year effect
  annual             <- dplyr::select(int, groupFctr, groupID, meanint, sdint)
  annual$smallPred   <- mainSmall+annual$meanint
  annual$largePred   <- mainLarge+annual$meanint
  annual$smallPredUp <- mainSmall+(annual$meanint+1.96*annual$sdint)
  annual$largePredUp <- mainLarge+(annual$meanint+1.96*annual$sdint)
  annual$smallPredDn <- mainSmall+(annual$meanint-1.96*annual$sdint)
  annual$largePredDn <- mainLarge+(annual$meanint-1.96*annual$sdint)
  annual$year        <- as.numeric(as.character(annual$groupID))

  ##  Subtract the cross-year average to calculate anomalies
  avg_small          <- mean(annual$smallPred)
  avg_large          <- mean(annual$largePred)
  annual$smallPred   <- annual$smallPred - avg_small
  annual$largePred   <- annual$largePred - avg_large
  annual$smallPredDn <- annual$smallPredDn - avg_small
  annual$smallPredUp <- annual$smallPredUp - avg_small
  annual$largePredDn <- annual$largePredDn - avg_large
  annual$largePredUp <- annual$largePredUp - avg_large


  ####
  ####  Anomalies for random year on intercept and slope model
  ####
  ##  Bootstrap effects
  fixed   <- FEsim(mod_rye_size_year, n.sims = 500) # bootstrapped fixed effects
  randoms <- REsim(mod_rye_size_year, n.sims = 500) # bootstrapped random effects

  ##  Calculate main effects for large and small plants
  mainSmall <- fixed$mean[1] + small*fixed$mean[2] + Wvals[1]*fixed$mean[3] + small*Wvals[1]*fixed$mean[4]
  mainLarge <- fixed$mean[1] + large*fixed$mean[2] + Wvals[2]*fixed$mean[3] + large*Wvals[2]*fixed$mean[4]

  ##  Extract random intercept parameters
  int <- randoms[randoms$term=="(Intercept)" & randoms$groupFctr=="year",]
  int$meanint <- int$mean
  int$sdint   <- int$sd
  slo <- randoms[randoms$term=="logarea.t0" & randoms$groupFctr=="year",]
  slo$meanslo <- slo$mean
  slo$sdslo   <-slo$sd
  annual_sxy <- merge(int[,-(3:6)],slo[,-(3:6)], by=c("groupID","groupFctr"))

  ##  Calculate the year effect
  annual_sxy$smallPred   <- mainSmall + annual_sxy$meanint + small*annual_sxy$meanslo
  annual_sxy$largePred   <- mainLarge + annual_sxy$meanint + large*annual_sxy$meanslo
  annual_sxy$smallPredUp <- mainSmall + annual_sxy$meanint + small*(annual_sxy$meanslo+1.96*annual_sxy$sdslo)
  annual_sxy$largePredUp <- mainLarge + annual_sxy$meanint + large*(annual_sxy$meanslo+1.96*annual_sxy$sdslo)
  annual_sxy$smallPredDn <- mainSmall + annual_sxy$meanint + small*(annual_sxy$meanslo-1.96*annual_sxy$sdslo)
  annual_sxy$largePredDn <- mainLarge + annual_sxy$meanint + large*(annual_sxy$meanslo-1.96*annual_sxy$sdslo)
  annual_sxy$year        <- as.numeric(as.character(annual_sxy$groupID))

  ##  Subtract the cross-year average to calculate anomalies
  avg_small_sxy          <- mean(annual_sxy$smallPred)
  avg_large_sxy          <- mean(annual_sxy$largePred)
  annual_sxy$smallPred   <- annual_sxy$smallPred - avg_small_sxy
  annual_sxy$largePred   <- annual_sxy$largePred - avg_large_sxy
  annual_sxy$smallPredDn <- annual_sxy$smallPredDn - avg_small_sxy
  annual_sxy$smallPredUp <- annual_sxy$smallPredUp - avg_small_sxy
  annual_sxy$largePredDn <- annual_sxy$largePredDn - avg_large_sxy
  annual_sxy$largePredUp <- annual_sxy$largePredUp - avg_large_sxy


  ####
  ####  Add ID information and store
  ####
  annual$state       <- do_state
  annual$species     <- do_species
  annual_sxy$state   <- do_state
  annual_sxy$species <- do_species

  rye_year_anoms <- rbind(rye_year_anoms, annual)
  rye_sxy_anoms  <- rbind(rye_sxy_anoms, annual_sxy)

  ##  Update progress
  cat(paste("Done with growth anomalies for",do_species,"in",do_state,"\n"))
  
} # end state-species growth loop

##  Save growth results
saveRDS(object = list(year_only = rye_year_anoms, size_by_year = rye_sxy_anoms),
        file = paste0(outDir,"grow_anomalies.RDS"))
write.csv(x = mod_comp_all,file = paste0(outDir, "grow_model_comparisons.csv"))

##  Remove stored objects before moving on
rm(rye_sxy_anoms,
   rye_year_anoms,
   mod_comp_all,
   mod_comp_tmp,
   gdata,
   gtmp,
   mod_rye_size_year,
   mod_rye_year,
   ll_full,
   ll_null,
   lrt,
   obs_lrt)



####
####  FIT SURVIVAL MODELS AND CALCULATE ANOMALIES --------------------------------
####
##  Define empty storage objects
rye_year_anoms <- {}
rye_sxy_anoms  <- {}
mod_comp_all   <- {}

## Loop over state-species
for(i in 1:nrow(state_species)){
  print(paste("Working on", as.character(state_species[i, 1]), as.character(state_species[i, 2])))
  ##  Set data directory
  do_state   <- state_species[i, "state"]
  do_species <- state_species[i, "species"]
  data_dir   <- paste0(datadrivs, do_state)

  ## Import demographic data
  dat <- fetchDemoData(do_species,
                       do_state,
                       grow=FALSE,
                       surv=TRUE,
                       recr=FALSE,
                       data_dir)
  sdata <- dat$data %>%
    dplyr::select(ID, year, logarea, survives, W, age, Group, quad) %>%
    mutate(year = as.factor(year))
  if(min(sdata$W) == 0){sdata$W <- sdata$W+.0001}

  ##  Drop seedlings and likely seedlings
  sdata <- trimQuadrats(data=sdata,
                        skip.years = TRUE,
                        df="df23",
                        trans="surv")$data

  ##  Fit models
  mod_rye_size_year <- glmer(survives~logarea*W+(1|Group)+(logarea.t0|year),
                             data=sdata,
                             family="binomial")
  mod_rye_year      <- glmer(survives~logarea*W+(1|Group)+(1|year),
                             data=sdata,
                             family="binomial")
  
  ## Loop over many iterations, simulate, compute LL: generates null distribution
  num_iters <- 500
  lrt <- numeric(num_iters)
  
  pb <- txtProgressBar(min = 1, max = num_iters, style = 3)
  for(iiter in 1:num_iters){
    stmp <- sdata # make copy of gdata
    y_star <- simulate(mod_rye_year) # simulate model
    stmp$survives <- as.numeric(y_star[,1])
    
    # Fit models with simulated data from the null hypothesis
    tmp_null <- glmer(survives~logarea.t0*W+(1|Group)+(1|year),
                      data = stmp, 
                      family = "binomial")
    tmp_full <- glmer(survives~logarea.t0*W+(1|Group)+(logarea.t0|year), 
                      data = stmp, 
                      family = "binomial")
    
    # Calculate log likelihoods
    ll_null <- logLik(tmp_null)
    ll_full <- logLik(tmp_full)
    lrt[iiter] <- 2*as.numeric(ll_full - ll_null)
    setTxtProgressBar(pb, iiter)
  }
  
  ## Calculate upper 95th quantile as value to exceed, for one-sided test
  upper <- as.numeric(quantile(lrt, probs = 0.95))
  
  ## Calculuate observed log likelihoods
  ll_null <- logLik(mod_rye_year)
  ll_full <- logLik(mod_rye_size_year)
  obs_lrt <- 2*as.numeric(ll_full - ll_null)
  
  ##  Store likelihood ratio tests
  mod_comp_tmp <- data.frame(state = do_state,
                             species = do_species,
                             obs_lrt = obs_lrt,
                             uquant = upper)
  mod_comp_all <- rbind(mod_comp_all, mod_comp_tmp)
  
  

  ####
  ####  Anomalies for random year on intercept model
  ####
  ##  Bootstrap effects
  fixed   <- FEsim(mod_rye_year, n.sims = 500) # bootstrapped fixed effects
  randoms <- REsim(mod_rye_year, n.sims = 500) # bootstrapped random effects

  ##  Calculate 10th and 90th percentiles of plant size
  small <- as.numeric(quantile(sdata$logarea.t0, probs=c(0.1)))
  large <- as.numeric(quantile(sdata$logarea.t0, probs=c(0.9)))

  ##  Calculate average crowding for large and small plants
  Wvals    <- c(NA,NA)
  Wvals[1] <- mean(sdata$W[which(sdata$logarea.t0 < small)])
  Wvals[2] <- mean(sdata$W[which(sdata$logarea.t0 > large)])

  ##  Calculate main effects for large and small plants
  mainSmall <- fixed$mean[1] + small*fixed$mean[2] + Wvals[1]*fixed$mean[3] + small*Wvals[1]*fixed$mean[4]
  mainLarge <- fixed$mean[1] + large*fixed$mean[2] + Wvals[2]*fixed$mean[3] + large*Wvals[2]*fixed$mean[4]

  ##  Extract random intercept parameters
  int <- randoms[randoms$term=="(Intercept)" & randoms$groupFctr=="year",]
  int$meanint <- int$mean
  int$sdint   <- int$sd

  ##  Calculate the year effect
  annual             <- dplyr::select(int, groupFctr, groupID, meanint, sdint)
  annual$smallPred   <- inv.logit(mainSmall+annual$meanint)
  annual$largePred   <- inv.logit(mainLarge+annual$meanint)
  annual$smallPredUp <- inv.logit(mainSmall+(annual$meanint+1.96*annual$sdint))
  annual$largePredUp <- inv.logit(mainLarge+(annual$meanint+1.96*annual$sdint))
  annual$smallPredDn <- inv.logit(mainSmall+(annual$meanint-1.96*annual$sdint))
  annual$largePredDn <- inv.logit(mainLarge+(annual$meanint-1.96*annual$sdint))
  annual$year        <- as.numeric(as.character(annual$groupID))

  ##  Subtract the cross-year average to calculate anomalies
  avg_small          <- mean(annual$smallPred)
  avg_large          <- mean(annual$largePred)
  annual$smallPred   <- annual$smallPred - avg_small
  annual$largePred   <- annual$largePred - avg_large
  annual$smallPredDn <- annual$smallPredDn - avg_small
  annual$smallPredUp <- annual$smallPredUp - avg_small
  annual$largePredDn <- annual$largePredDn - avg_large
  annual$largePredUp <- annual$largePredUp - avg_large


  ####
  ####  Anomalies for random year on intercept and slope model
  ####
  ##  Bootstrap effects
  fixed   <- FEsim(mod_rye_size_year, n.sims = 500) # bootstrapped fixed effects
  randoms <- REsim(mod_rye_size_year, n.sims = 500) # bootstrapped random effects

  ##  Calculate main effects for large and small plants
  mainSmall <- fixed$mean[1] + small*fixed$mean[2] + Wvals[1]*fixed$mean[3] + small*Wvals[1]*fixed$mean[4]
  mainLarge <- fixed$mean[1] + large*fixed$mean[2] + Wvals[2]*fixed$mean[3] + large*Wvals[2]*fixed$mean[4]

  ##  Extract random intercept parameters
  int <- randoms[randoms$term=="(Intercept)" & randoms$groupFctr=="year",]
  int$meanint <- int$mean
  int$sdint   <- int$sd
  slo <- randoms[randoms$term=="logarea.t0" & randoms$groupFctr=="year",]
  slo$meanslo <- slo$mean
  slo$sdslo   <-slo$sd
  annual_sxy <- merge(int[,-(3:6)],slo[,-(3:6)], by=c("groupID","groupFctr"))

  ##  Calculate the year effect
  annual_sxy$smallPred   <- inv.logit(mainSmall + annual_sxy$meanint + small*annual_sxy$meanslo)
  annual_sxy$largePred   <- inv.logit(mainLarge + annual_sxy$meanint + large*annual_sxy$meanslo)
  annual_sxy$smallPredUp <- inv.logit(mainSmall + annual_sxy$meanint + small*(annual_sxy$meanslo+1.96*annual_sxy$sdslo))
  annual_sxy$largePredUp <- inv.logit(mainLarge + annual_sxy$meanint + large*(annual_sxy$meanslo+1.96*annual_sxy$sdslo))
  annual_sxy$smallPredDn <- inv.logit(mainSmall + annual_sxy$meanint + small*(annual_sxy$meanslo-1.96*annual_sxy$sdslo))
  annual_sxy$largePredDn <- inv.logit(mainLarge + annual_sxy$meanint + large*(annual_sxy$meanslo-1.96*annual_sxy$sdslo))
  annual_sxy$year        <- as.numeric(as.character(annual_sxy$groupID))

  ##  Subtract the cross-year average to calculate anomalies
  avg_small_sxy          <- mean(annual_sxy$smallPred)
  avg_large_sxy          <- mean(annual_sxy$largePred)
  annual_sxy$smallPred   <- annual_sxy$smallPred - avg_small_sxy
  annual_sxy$largePred   <- annual_sxy$largePred - avg_large_sxy
  annual_sxy$smallPredDn <- annual_sxy$smallPredDn - avg_small_sxy
  annual_sxy$smallPredUp <- annual_sxy$smallPredUp - avg_small_sxy
  annual_sxy$largePredDn <- annual_sxy$largePredDn - avg_large_sxy
  annual_sxy$largePredUp <- annual_sxy$largePredUp - avg_large_sxy


  ####
  ####  Add ID information and store
  ####
  annual$state       <- do_state
  annual$species     <- do_species
  annual_sxy$state   <- do_state
  annual_sxy$species <- do_species
  mod_comp$state     <- do_state
  mod_comp$species   <- do_species

  rye_year_anoms <- rbind(rye_year_anoms, annual)
  rye_sxy_anoms  <- rbind(rye_sxy_anoms, annual_sxy)

  ##  Update progress
  cat(paste("Done with survival anomalies for",do_species,"in",do_state,"\n"))

} # end state-species growth loop

##  Save survival results
saveRDS(object = list(year_only = rye_year_anoms, size_by_year = rye_sxy_anoms),
        file = paste0(outDir,"surv_anomalies.RDS"))
write.csv(x = mod_comp_all,file = paste0(outDir, "surv_model_comparisons.csv"))


