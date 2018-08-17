################################################################################
##  plot_expected_dists.R: Script to plot the distribution of expected sizes
##  given current size and model type (size-by-year effect, or no).
##
##  Author: Andrew Tredennick (atredenn@gmail.com)
##  Date created: April 16, 2018
################################################################################

rm(list = ls(all.names = TRUE)) # clear the workspace



####
####  LOAD LIBRARIES -----------------------------------------------------------
####
library(tidyverse) # Data manipulation
library(stringr)   # Workging with strings
library(ggthemes)  # Pleasing ggplot2 themes
library(ggalt)     # Some more geoms for ggplot2
library(cowplot)   # Combining ggplot2 objects



####
####  SET DIRECTORIES ----------------------------------------------------------
####
root        <- ifelse(.Platform$OS.type=="windows","c:/repos","~/repos")
work_dir    <- paste0(root,"/drivers/empirical/size_by_year_IPM/analysis/")
results_dir <- paste0(root,"/drivers/emprical/size_by_year_IPM/")
ms_fig_dir <- paste0(root,"/drivers/Manuscripts/SizebyClimate/figures/")
setwd(work_dir)



####
####  DEFINE USER FUNCTIONS ----------------------------------------------------
####
##  Calculate expectations from regressions
get_expectations <- function(dofile, filepath, model_to_use, zpreds){
  tmp_model <- model_to_use
  tmp_vital <- str_split(dofile, pattern = "_")[[1]][1]
  tmp_species <- str_split(dofile, pattern = "_")[[1]][4]
  tmp_species <- str_split(tmp_species, pattern = "[.]")[[1]][1]
  tmp_state <- str_split(dofile, pattern = "_")[[1]][3]
  tmp_params <- readRDS(paste0(filepath,dofile))[[model_to_use]]
  
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
  out_df$model <- model_to_use
  return(out_df)
}

##  Function for getting empirical plant size distribution
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

##  Antilogit
antilogit <- function(x) { exp(x) / (1 + exp(x) ) }


####
####  GRAB GROWTH REGRESSION PARAMETERS AND PREDICT ----------------------------
####
##  Settings
zpreds <- seq(-2,10,by=0.1) # sizes at which to make predictions
model_types <- c("yearonly","sizeXyear")

##  One-stage IPM settings
onestage_files <- list.files(path = "../one_stage/vital_rate_params/")
recones <- grep("rec", onestage_files)
onestage_files <- onestage_files[-recones]
onestagers <- data.frame(path = "../one_stage/vital_rate_params/",
                         files = onestage_files)

##  Two-stage IPM settings
twostage_files <- list.files(path = "../two_stage/vital_rate_params/")
recones <- grep("rec", twostage_files)
twostage_files <- twostage_files[-recones]
twostagers <- data.frame(path = "../two_stage/vital_rate_params/",
                         files = twostage_files)

allfiles <- rbind(rbind(onestagers, twostagers),
                  rbind(onestagers, twostagers))
allfiles$modeltype <- rep(model_types, each = (nrow(allfiles)/2))



####
####  LOOP THROUGH FILES AND MODEL TYPES TO PREDICT SIZES ----------------------
####
small_large_preds <- {}
for(i in 1:nrow(allfiles)){
  ipath <- allfiles$path[i]
  ifile <- allfiles$files[i]
  imodel <- allfiles$modeltype[i]
  
  preds <- get_expectations(ifile, ipath, imodel, zpreds)
  
  plant_sizes <- get_plant_sizes(state = unique(preds$state), 
                                 species = unique(preds$species))
  
  # Get the closest real values to the quantile estimates for large and small
  zL_real <- preds$z0[which.min(abs(preds$z0 - plant_sizes[1]))]
  zU_real <- preds$z0[which.min(abs(preds$z0 - plant_sizes[2]))]
  
  small_large_sub <- preds %>%
    filter(z0 %in% c(zL_real, zU_real)) %>%
    mutate(size = ifelse(z0 == zL_real, "Small plant", "Large plant"))
  small_large_preds <- rbind(small_large_preds, small_large_sub)
}



####
####  CALCULATE AVERAGES IN LOG AND REAL SPACE ---------------------------------
####
growth_means <- small_large_preds %>%
  filter(vital == "growth") %>%
  mutate(exp_z1_hat = exp(z1_hat)) %>%
  group_by(state, species, vital, size, model) %>%
  summarise(mean_z1 = mean(z1_hat),
            mean_expz1 = mean(exp_z1_hat)) %>%
  gather(transformation, avg_pred, mean_z1:mean_expz1) %>%
  spread(model, avg_pred) %>%
  mutate(trans_name = ifelse(transformation == "mean_z1", "Log scale", "Arithmetic scale"))

surv_means <- small_large_preds %>%
  filter(vital == "surv") %>%
  mutate(prob_z1_hat = antilogit(z1_hat)) %>%
  group_by(state, species, vital, size, model) %>%
  summarise(mean_z1 = mean(z1_hat),
            mean_probz1 = mean(prob_z1_hat)) %>%
  gather(transformation, avg_pred, mean_z1:mean_probz1) %>%
  spread(model, avg_pred) %>%
  mutate(trans_name = ifelse(transformation == "mean_z1", "Logit scale", "Antilogit scale"))

  
ggplot(growth_means, aes(x = yearonly, y = sizeXyear))+
  geom_abline(aes(intercept = 0, slope = 1), color = "grey")+
  geom_point(size = 2)+
  facet_wrap(size~trans_name, scales = "free",
             labeller = function (labels) {
               labels <- lapply(labels, as.character)
               list(do.call(paste, c(labels, list(sep = "\n"))))
             })+
  xlab("Mean E(z) from size-independent model")+
  ylab("Mean E(z) from size-dependent model")+
  theme_few()
ggsave("../../../Manuscripts/SizebyClimate/figures/grow_expectations.pdf", 
       width = 6, 
       height = 6, 
       units = "in")

ggplot(surv_means, aes(x = yearonly, y = sizeXyear))+
  geom_abline(aes(intercept = 0, slope = 1), color = "grey")+
  geom_point(size = 2)+
  facet_wrap(size~trans_name, scales = "free",
             labeller = function (labels) {
               labels <- lapply(labels, as.character)
               list(do.call(paste, c(labels, list(sep = "\n"))))
             })+
  xlab("Mean Pr(S) from size-independent model")+
  ylab("Mean Pr(S) from size-dependent model")+
  theme_few()
ggsave("../../../Manuscripts/SizebyClimate/figures/surv_expectations.pdf", 
       width = 6, 
       height = 6, 
       units = "in")

# 
# ####
# ####  MAKE THE PLOTS -----------------------------------------------------------
# ####
# growth <- small_large_preds %>%
#   filter(vital == "growth") %>%
#   mutate(state_spp = paste(state, species, sep = ", "),
#          exp_z1_hat = exp(z1_hat)) %>%
#   select(year, state_spp, model, size, z1_hat, exp_z1_hat) %>%
#   gather(key = "transformation", value = "value", z1_hat:exp_z1_hat)
# 
# survival <- small_large_preds %>%
#   filter(vital == "surv") %>%
#   mutate(state_spp = paste(state, species, sep = ", "),
#          prob_z1_hat = antilogit(z1_hat)) %>%
#   select(year, state_spp, model, size, z1_hat, prob_z1_hat) %>%
#   gather(key = "transformation", value = "value", z1_hat:prob_z1_hat)
# 
# # Function for making and saving the plots
# make_save_plot <- function(df, fname, plant_size, transform_str){
#   logplot <- ggplot(filter(df, size == plant_size & transformation == transform_str[1]),
#                       aes(x = value))+
#     geom_density(aes(color = state_spp, linetype = model))+
#     facet_wrap(~state_spp, scales = "free", ncol = 5)+
#     guides(color = FALSE)+
#     theme_grey()+
#     theme(axis.text.y = element_blank())+
#     ggtitle("Log scale expectations")
#   
#   arithplot <- ggplot(filter(df, size == plant_size & transformation == transform_str[2]),
#                       aes(x = value))+
#     geom_density(aes(color = state_spp, linetype = model))+
#     facet_wrap(~state_spp, scales = "free", ncol = 5)+
#     guides(color = FALSE)+
#     theme_grey()+
#     theme(axis.text.y = element_blank())+
#     ggtitle("Natural scale expectations")
#   
#   theplots <- plot_grid(logplot, arithplot, 
#                             nrow = 2,
#                             align = "v",
#                             labels = "AUTO")
#   ggsave(filename = fname, 
#          plot = theplots, height = 11, width = 8, units = "in")
# }
# 
# make_save_plot(df = growth, 
#                fname = "/Users/atredenn/Desktop/growth_large.png", 
#                plant_size = "Large plant", 
#                transform_str = unique(growth$transformation))
# 
# make_save_plot(df = growth, 
#                fname = "/Users/atredenn/Desktop/growth_small.png", 
#                plant_size = "Small plant", 
#                transform_str = unique(growth$transformation))
# 
# make_save_plot(df = survival, 
#                fname = "/Users/atredenn/Desktop/surv_large.png", 
#                plant_size = "Large plant", 
#                transform_str = unique(survival$transformation))
# 
# make_save_plot(df = survival, 
#                fname = "/Users/atredenn/Desktop/surv_small.png", 
#                plant_size = "Small plant", 
#                transform_str = unique(survival$transformation))
# 
