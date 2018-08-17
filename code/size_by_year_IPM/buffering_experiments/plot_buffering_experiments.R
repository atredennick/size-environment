################################################################################
##  plot_buffering_experiments.R: Script to plot the distributions of plant
##  cover from the buffering experiment simulations.
##
##  Author: Andrew Tredennick (atredenn@gmail.com)
##  Date created: 11-6-2017
################################################################################

rm(list = ls(all.names = TRUE))

####
####  LOAD LIBRARIES ----
####
library(tidyverse)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggthemes)
library(viridis)



####
####  SET DIRECTORIES ----
####
root     <- ifelse(.Platform$OS.type=="windows","c:/repos","~/repos")
work_dir <- paste0(root,"/drivers/empirical/size_by_year_IPM/buffering_experiments/")
data_dir <- paste0(root,"/driversdata/data/")
res_dir  <- paste0(root,"/drivers/empirical/size_by_year_IPM/buffering_experiments/results/")
ms_fig_dir  <- paste0(root, "/drivers/Manuscripts/SizebyClimate/figures/")
setwd(work_dir) # set working directory



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

all_sims <- {}
for(doit in 1:nrow(state_species)){
  do_state      <- as.character(state_species[doit, "state"])
  do_species <- as.character(state_species[doit, "species"])
  
  transform_matrix <- expand.grid(c("variable_small","variable_large","sign_flip"),
                                  c("variable_small","variable_large","sign_flip"))
  transform_matrix <- rbind(transform_matrix,
                            data.frame(Var1 = c("variable_small","variable_large","sign_flip"),
                                       Var2 = rep("none",3)),
                            data.frame(Var1 = rep("none",3),
                                       Var2 = c("variable_small","variable_large","sign_flip")),
                            data.frame(Var1 = "none", Var2 = "none"))
  colnames(transform_matrix) <- c("growth","survival")
  transform_matrix$sim_id <- as.character(1:nrow(transform_matrix))
  
  res_dir  <- paste0(root,"/drivers/empirical/size_by_year_IPM/buffering_experiments/lambda_results/")
  
  all_files <- as.data.frame(list.files(res_dir)) %>%
    rename(fname = `list.files(res_dir)`) %>%
    separate(fname, into = c("state", "species","simext"), sep = "_") %>%
    separate(simext, into = c("sim_id","ext"), sep = "[.]") %>%
    cbind(as.data.frame(list.files(res_dir))) %>%
    rename(fname = `list.files(res_dir)`) %>%
    dplyr::select(fname, state, species, sim_id) %>%
    left_join(transform_matrix, by = "sim_id") %>%
    mutate(sim_id = as.numeric(sim_id)) %>%
    arrange(sim_id) %>%
    filter(state == do_state & species == do_species)
  
  # all_sims <- {} # create empty object for storage
  for(i in 1:length(all_files$fname)){
    dofile <- all_files$fname[i]
    tmp <- readRDS(paste0("./lambda_results/",dofile))
    if("cov_kids" %in% colnames(tmp)){
      tmp <- tmp %>%
        dplyr::mutate(cover = cov_kids + cov_adults) %>%
        dplyr::select(iteration, year, cover)
    }
    start_cover <- tmp$cover[1]
    growth_rates <- log(tmp$cover[2:nrow(tmp)] / start_cover)
    mean_pgr <- mean(growth_rates)
    tmp_out <- data.frame(sim_id = all_files$sim_id[i],
                          lambda = mean_pgr,
                          species = do_species,
                          state = do_state)
    all_sims <- rbind(all_sims, tmp_out)
  } # end simulation loop
} # end state, species loop


sims_summ <- all_sims %>%
  dplyr::group_by(state, species) %>%
  dplyr::mutate(null_lambda = lambda[sim_id == 16], # null lamda is the no-perturbation simulation
                lambda_scaled = lambda / null_lambda) %>%
  dplyr::ungroup() %>%
  mutate(sim_id = as.character(sim_id)) %>%
  left_join(transform_matrix, by = "sim_id") %>%
  mutate(simulation = paste0("Growth: ", growth,"\nSurvival: ", survival),
         colorid = ifelse(sim_id == 16, "a", "b"),
         state_spp = paste(state, species, sep = ", "))


####
####  MAKE THE PLOT ------------------------------------------------------------
####
better_names <- data.frame(state = c("idaho","montana","arizona","newmexico","kansas"),
                           new_state = c("Idaho","Montana","Arizona","New Mexico","Kansas")) %>%
  dplyr::mutate(state = as.character(state),
                new_state = as.character(new_state))

plot_summ <- sims_summ %>%
  left_join(better_names, by = c("state")) %>%
  dplyr::mutate(state_spp = paste(new_state, species, sep = ", ")) %>%
  filter(growth != "sign_flip") %>%
  filter(survival != "sign_flip")



mylabs <- c(
  expression(paste("Low ", lambda[S])),
  "",
  "",
  expression(paste(lambda[S], " = 1")),
  "",
  "",
  expression(paste("High ", lambda[S]))
)

buffs <- ggplot(plot_summ, aes(x=growth, y=survival, fill = lambda_scaled))+
  geom_tile()+
  geom_text(aes(label = round(lambda_scaled, 2)), color = "grey25")+
  scale_fill_viridis(name = expression(paste("Standardized ", lambda[S])), limits = c(0.7,1.31), labels = mylabs)+
  scale_x_discrete(labels = c("Small more variable","Large more variable","None"))+
  scale_y_discrete(labels = c("Small more variable","Large more variable","None"))+
  ylab("Survival Effect")+
  xlab("Growth Effect")+
  coord_equal()+
  facet_wrap(~state_spp, ncol = 5)+
  theme_few()+
  theme(axis.text.x = element_text(angle = 45, hjust =1))
ggsave(paste0(ms_fig_dir,"buffering_results.pdf"), plot = buffs, width = 9.5, height = 6, units = "in")

