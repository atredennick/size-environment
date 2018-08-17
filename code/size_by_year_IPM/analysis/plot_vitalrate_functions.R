################################################################################
##  plot_vitalrate_functions.R: Script that plots the empirical vital rate
##  functions for survival and growth on log and aritmetics scales.
##
##  Author: Andrew Tredennick (atredenn@gmail.com)
##  Date created: April 11, 2018
################################################################################

rm(list = ls(all.names = TRUE)) # clear the workspace



####
####  LOAD LIBRARIES -----------------------------------------------------------
####
library(tidyverse) # Data manipulation
library(stringr)   # Workging with strings
library(ggthemes)  # Pleasing ggplot2 themes
library(ggalt)     # Some more geoms for ggplot2
# library(cowplot)   # Combining ggplot2 objects



####
####  SET DIRECTORIES ----------------------------------------------------------
####
root        <- ifelse(.Platform$OS.type=="windows","c:/repos","~/repos")
work_dir    <- paste0(root,"/drivers/empirical/size_by_year_IPM/analysis/")
results_dir <- paste0(root,"/drivers/emprical/size_by_year_IPM/")
ms_fig_dir <- paste0(root,"/drivers/Manuscripts/SizebyClimate/figures/")
setwd(work_dir)


####
####  GRAB GROWTH REGRESSION PARAMETERS AND PREDICT ----------------------------
####
##  Settings
zpreds <- seq(-2,10,by=0.1) # sizes at which to make predictions
model_to_use <- "yearonly"
# model_to_use <- "sizeXyear"

##  One-stage IPM
onestage_files <- list.files(path = "../one_stage/vital_rate_params/")
recones <- grep("rec", onestage_files)
onestage_files <- onestage_files[-recones]

allpreds_df <- {}
for(dofile in onestage_files){
  tmp_vital <- str_split(dofile, pattern = "_")[[1]][1]
  tmp_species <- str_split(dofile, pattern = "_")[[1]][4]
  tmp_species <- str_split(tmp_species, pattern = "[.]")[[1]][1]
  tmp_state <- str_split(dofile, pattern = "_")[[1]][3]
  tmp_params <- readRDS(paste0("../one_stage/vital_rate_params/",dofile))[[model_to_use]]

  
  year_id <- expand.grid(tmp_params$year, zpreds)
  intercept_eff <- expand.grid(tmp_params$intercept, zpreds)
  if(tmp_vital == "growth"){
    slope_eff <- expand.grid(tmp_params$logarea.t0, zpreds)
  }
  if(tmp_vital == "surv"){
    slope_eff <- expand.grid(tmp_params$logarea, zpreds)
  }
  zouts <- intercept_eff[,1] + slope_eff[,1]*slope_eff[,2]
 
  out_df <- cbind(year_id, zouts)
  colnames(out_df) <- c("year","z0","z1_hat")
  out_df$species <- tmp_species
  out_df$state <- tmp_state
  out_df$vital <- tmp_vital
  
  allpreds_df <- rbind(allpreds_df, out_df)
}

##  Two-stage IPM
twostage_files <- list.files(path = "../two_stage/vital_rate_params/")
recones <- grep("rec", twostage_files)
twostage_files <- twostage_files[-recones]

for(dofile in twostage_files){
  tmp_vital <- str_split(dofile, pattern = "_")[[1]][1]
  tmp_species <- str_split(dofile, pattern = "_")[[1]][4]
  tmp_species <- str_split(tmp_species, pattern = "[.]")[[1]][1]
  tmp_state <- str_split(dofile, pattern = "_")[[1]][3]
  tmp_params <- readRDS(paste0("../two_stage/vital_rate_params/",dofile))[[model_to_use]]
  
  year_id <- expand.grid(tmp_params$year, zpreds)
  intercept_eff <- expand.grid(tmp_params$intercept, zpreds)
  if(tmp_vital == "growth"){
    slope_eff <- expand.grid(tmp_params$logarea.t0, zpreds)
  }
  if(tmp_vital == "surv"){
    slope_eff <- expand.grid(tmp_params$logarea, zpreds)
  }
  zouts <- intercept_eff[,1] + slope_eff[,1]*slope_eff[,2]
  
  out_df <- cbind(year_id, zouts)
  colnames(out_df) <- c("year","z0","z1_hat")
  out_df$species <- tmp_species
  out_df$state <- tmp_state
  out_df$vital <- tmp_vital
  
  allpreds_df <- rbind(allpreds_df, out_df)
}



####
####  PLOT EMPIRICAL GROWTH AND SURVIVAL FUNCTIONS -----------------------------
####
growth_preds <- allpreds_df %>%
  filter(vital == "growth") %>%
  mutate(state_spp = paste(state, species, sep = ", "))

antilogit <- function(x) { exp(x) / (1 + exp(x) ) }
surv_preds <- allpreds_df %>%
  filter(vital == "surv") %>%
  mutate(state_spp = paste(state, species, sep = ", ")) %>%
  mutate(sprob = antilogit(z1_hat))

ggplot(growth_preds, aes(x = exp(z0), y = exp(z1_hat), color = as.factor(year)))+
  geom_line()+
  guides(color = FALSE)+
  facet_wrap(~state_spp, scales = "free")+
  theme_few()+
  xlab(expression(paste("Observed size, ",z[t])))+
  ylab(expression(paste("Predicted size, ",z[t+1])))

ggplot(surv_preds, aes(x = exp(z0), y = sprob, color = as.factor(year)))+
  geom_line()+
  guides(color = FALSE)+
  facet_wrap(~state_spp, scales = "free")+
  theme_few()+
  xlab(expression(paste("Observed size, ",z[t])))+
  ylab("Predicted survival, Pr(S)")



####
####  GRAB SLICES FROM SMALL AND LARGE PLANTS AND PLOT HISTOGRAMS --------------
####
get_plant_sizes <- function(state, species){
  state <- tolower(state)
  species <- toupper(species)
  basedir <- "../../../../driversdata/data/"
  all_files <- list.files(paste0(basedir,state,"/speciesData/",species,"/"))
  survfiles <- all_files[grep("surv", all_files)]
  survfile <- survfiles[length(survfiles)]
  obs_data <- read.csv(paste0(basedir,state,"/speciesData/",species,"/",survfile))
  logareas <- log(obs_data$area)
  small <- quantile(logareas, probs = 0.05)
  large <- quantile(logareas, probs = 0.95)
  return(c(small,large))
} 

dostate <- "kansas"
dospecies <- "ANGE"
plant_sizes <- get_plant_sizes(state = dostate, species = dospecies)
tmp_preds <- filter(growth_preds, state == dostate & species == dospecies)

# Get the closest real values to the quantile estimates for large and small
zL_real <- tmp_preds$z0[which.min(abs(tmp_preds$z0 - plant_sizes[1]))]
zU_real <- tmp_preds$z0[which.min(abs(tmp_preds$z0 - plant_sizes[2]))]

# Subset out only large and small plants
small_large_sub <- tmp_preds %>%
  filter(z0 %in% c(zL_real, zU_real)) %>%
  mutate(size = ifelse(z0 == zL_real, "Small plant", "Large plant"))

# Plot their yhat distributions
ggplot(small_large_sub, aes(x = exp(z1_hat), fill = size))+
  geom_density(color = NA)+
  facet_wrap(~size, scales = "free")+
  scale_fill_brewer(type = "qual")+
  theme_few()+
  guides(fill = FALSE)+
  xlab(expression(paste("Predicted size, ",z[t+1])))+
  ylab("Density")+
  ggtitle("Distribution of predicted sizes: ANGE (Kansas)")

## Now do all species
all_small_large <- {}
for(dostate in unique(growth_preds$state)){
  tmpdf <- filter(growth_preds, state == dostate)
  for(dospecies in unique(tmpdf$species)){
    plant_sizes <- get_plant_sizes(state = dostate, species = dospecies)
    tmp_preds <- filter(growth_preds, state == dostate & species == dospecies)
    
    # Get the closest real values to the quantile estimates for large and small
    zL_real <- tmp_preds$z0[which.min(abs(tmp_preds$z0 - plant_sizes[1]))]
    zU_real <- tmp_preds$z0[which.min(abs(tmp_preds$z0 - plant_sizes[2]))]
  
    # Subset out only large and small plants
    small_large_sub <- tmp_preds %>%
      filter(z0 %in% c(zL_real, zU_real)) %>%
      mutate(size = ifelse(z0 == zL_real, "Small plant", "Large plant"))
    all_small_large <- rbind(all_small_large, small_large_sub)
  }
}

ggplot(all_small_large, aes(x = exp(z1_hat), fill = size))+
  geom_density(color = NA)+
  facet_wrap(state_spp~size, scales = "free")+
  scale_fill_brewer(type = "qual")+
  theme_few()+
  guides(fill = FALSE)+
  xlab(expression(paste("Predicted size, ",z[t+1])))+
  ylab("Density")+
  theme(axis.text.y = element_blank())



####
####  PLOT RANDOM INTERCEPT DISTRIBUTIONS --------------------------------------
####
model_to_use <- "yearonly"

##  One-stage IPM
onestage_files <- list.files(path = "../one_stage/vital_rate_params/")
recones <- grep("rec", onestage_files)
onestage_files <- onestage_files[-recones]

allpreds_df <- {}
for(dofile in onestage_files){
  tmp_vital <- str_split(dofile, pattern = "_")[[1]][1]
  tmp_species <- str_split(dofile, pattern = "_")[[1]][4]
  tmp_species <- str_split(tmp_species, pattern = "[.]")[[1]][1]
  tmp_state <- str_split(dofile, pattern = "_")[[1]][3]
  tmp_params <- readRDS(paste0("../one_stage/vital_rate_params/",dofile))[[model_to_use]]
  
  out_df <- data.frame(year = tmp_params$year_id,
                       intercept = tmp_params$intercept,
                       state = tmp_state,
                       species = tmp_species,
                       vital = tmp_vital)
  
  
  allpreds_df <- rbind(allpreds_df, out_df)
}

##  Two-stage IPM
twostage_files <- list.files(path = "../two_stage/vital_rate_params/")
recones <- grep("rec", twostage_files)
twostage_files <- twostage_files[-recones]

for(dofile in twostage_files){
  tmp_vital <- str_split(dofile, pattern = "_")[[1]][1]
  tmp_species <- str_split(dofile, pattern = "_")[[1]][4]
  tmp_species <- str_split(tmp_species, pattern = "[.]")[[1]][1]
  tmp_state <- str_split(dofile, pattern = "_")[[1]][3]
  tmp_params <- readRDS(paste0("../two_stage/vital_rate_params/",dofile))[[model_to_use]]
  
  out_df <- data.frame(year = tmp_params$year_id,
                       intercept = tmp_params$intercept,
                       state = tmp_state,
                       species = tmp_species,
                       vital = tmp_vital)
  
  allpreds_df <- rbind(allpreds_df, out_df)
}

growth_intercept <- allpreds_df %>%
  filter(vital == "growth") %>%
  mutate(state_spp = paste(state, species, sep = ", "),
         expint = exp(intercept)) %>%
  select(year, state_spp, intercept, expint) %>%
  gather(key = "transformation", value = "value", intercept:expint)

ggplot(growth_intercept, aes(x=value, fill = transformation))+
  geom_density(color = NA, alpha = 0.8)+
  facet_wrap(~state_spp, scales = "free")+
  scale_fill_brewer(type = "qual", name = NULL, labels = c("exp(intercept)", "intercept"))+
  theme_few()

# intercepts <- rnorm(1000,0.6,1)
# par(mfrow = c(2,1))
# hist(intercepts)
# abline(v = 0.6, col = "blue", lwd = 2)
# hist(exp(intercepts))
# abline(v = exp(0.6), col = "blue", lwd = 2)





