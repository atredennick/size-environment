################################################################################
##  analyze_ipm_results.R: script to plot results from IPM simulations with
##  and without size-by-year effects. The script produces boxplots of
##  equilibrium runs, time series plots of transient runs, and a table summary
##  of average time to reach equilibrium cover with and without the size*year
##  interaction.
##
##  Author: Andrew Tredennick (atredenn@gmail.com)
##  Date created: 10-11-2017
################################################################################

##  Clear everything
rm(list = ls(all.names = TRUE))


####
####  LOAD LIBRARIES ----
####
library(tidyverse) # Data manipulation
library(stringr)   # Workging with strings
library(ggthemes)  # Pleasing ggplot2 themes
library(ggalt)     # Some more geoms for ggplot2
library(cowplot)   # Combining ggplot2 objects



####
####  SET DIRECTORIES ----
####
root        <- ifelse(.Platform$OS.type=="windows","c:/repos","~/repos")
work_dir    <- paste0(root,"/drivers/empirical/size_by_year_IPM/analysis/")
results_dir <- paste0(root,"/drivers/empirical/size_by_year_IPM/ipm_results/")
loc_fig_dir <- paste0(root,"/drivers/empirical/size_by_year_IPM/analysis/figures/")
ms_fig_dir <- paste0(root,"/drivers/Manuscripts/SizebyClimate/figures/")


setwd(work_dir) # set working directory
dir.create(file.path(loc_fig_dir), showWarnings = FALSE) # create a folder for figures
dir.create(file.path(ms_fig_dir), showWarnings = FALSE) # create a folder for figures



####
####  READ IN AND COMBINE RESULTS FOR EACH SPECIES AND STATE ----
####
all_files <- list.files(paste0(results_dir,"equilibrium_runs/"))
all_sims  <- {}
for(do_file in all_files){
  file_info <- unlist(str_split(do_file,"_"))
  state     <- file_info[1]
  species   <- file_info[2]
  simtype   <- unlist(str_split(file_info[3],"[.]"))[1]
  
  ##  Rename states with formatting for figures
  if(state == "arizona"){state <- "Arizona"}
  if(state == "kansas"){state <- "Kansas"}
  if(state == "montana"){state <- "Montana"}
  if(state == "newmexico"){state <- "New Mexico"}
  if(state == "idaho"){state <- "Idaho"}
  
  tmp <- readRDS(paste0(results_dir,"equilibrium_runs/",do_file)) %>%
    mutate(state = state, species = species, simtype = simtype)
  
  if("cov_kids" %in% colnames(tmp)){
    tmp <- tmp %>%
      mutate(cover = cov_kids + cov_adults,
             abundance = N_kids + N_adults) %>%
      dplyr::select(iteration,year,cover,abundance,state,species,simtype)
    }
  
  all_sims <- rbind(all_sims, tmp)
}

cut_sims <- filter(all_sims, iteration>199)



####
####  CALCULATE STD. DEV. OF EACH RUN ----
####
summary_sims <- cut_sims %>%
  group_by(state, species, simtype) %>%
  summarise(sd_cover = sd(cover),
            mean_cover = mean(cover)) %>%
  mutate(sim_text = ifelse(simtype=="sizeXyear", "with SxY", "without SxY")) %>%
  mutate(sd_text = paste("S.D.", sim_text, "=",round(sd_cover,3))) %>%
  ungroup()
saveRDS(summary_sims, "../analysis/ipm_summary_stats.RDS")



####
####  MAKE THE PLOTS ----
####
equilibrium <- ggplot(cut_sims, aes(x=species,y=cover*100,fill=simtype))+
  geom_boxplot(width=0.5, outlier.size = 0.1, lwd=0.2)+
  #scale_fill_brewer(palette = "Set1", name = "", labels=c("Size*Year", "No Size*Year"))+
  # scale_fill_manual(values=c("grey50","white"), name = "", labels=c("Size*Year", "No Size*Year"))+
  scale_fill_brewer(type = "qual", labels=c("Size*Year", "No Size*Year"), name = "")+
  ylab("Cover (%)")+
  xlab("Species")+
  facet_wrap(~state, scales = "free", nrow=1)+
  theme_few()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# ggsave(paste0(loc_fig_dir,"ipm_boxplots.png"), plot = equilibrium, width = 9, height = 2.5, units = "in")

##  Time series
# ggplot(filter(cut_sims, iteration<251), aes(x=iteration, y = cover*100, color=species))+
#   geom_line(aes(linetype=simtype))+
#   facet_wrap(~state, scales = "free")
# # ggsave(paste0(out_dir,"ipm_timeseries.pdf"), width = 9, height = 5, units = "in")
# 
# ggplot(filter(cut_sims, iteration<251 & state=="Montana"), aes(x=iteration, y = cover*100, color=species))+
#   geom_line(aes(linetype=simtype))+
#   facet_wrap(~species, scales = "free")
# 
# ggplot(filter(cut_sims, iteration<251 & state=="Arizona"), aes(x=iteration, y = cover*100, color=species))+
#   geom_line(aes(linetype=simtype))+
#   facet_wrap(~species, scales = "free")



####
####  LOOK AT VARIANCE OF RANDOM YEAR EFFECTS -----
####
ipm_calls <- read.table("../run_all_ipms.sh") %>%
  dplyr::select(-V1) %>%
  dplyr::rename(path=V2, id=V3) %>%
  separate(path,c("tmp1","folder","tmp2","tmp3"),"/") %>%
  dplyr::select(folder,id)

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
state_species$size_by_year <- "sizeXyear"
state_species <- rbind(state_species,
                       mutate(state_species, size_by_year="yearonly"))

state_species <- cbind(state_species,ipm_calls)

states <- unique(state_species$state)
vitals <- c("growth","surv")
year_eff_df <- {} # empty object
for(do_state in states){
  species <- unique(filter(state_species, state==do_state)$species)
  for(do_species in species){
    folder <- unique(filter(state_species, state==do_state & species==do_species)$folder)
    for(do_vital in vitals){
      file_name <- paste0(do_vital,"_params_",do_state,"_",do_species,".RDS")
      tmp_params <- readRDS(paste0("../",folder,"/vital_rate_params/",file_name))[["sizeXyear"]]
      intercept_sd <- sd(tmp_params$intercept)
      if(do_vital == "growth"){
        size_sd <- sd(tmp_params$logarea.t0)
        effect_cor <- cor(tmp_params$intercept, tmp_params$logarea.t0)
      }
      if(do_vital == "surv"){
        size_sd <- sd(tmp_params$logarea)
        effect_cor <- cor(tmp_params$intercept, tmp_params$logarea)
      }
      # plot(tmp_params$intercept, tmp_params$logarea.t0, xlab="Intercept", ylab="logarea.t0",main=do_species)
      tmp_df <- data.frame(state = do_state,
                           species = do_species,
                           vital_rate = do_vital,
                           sd_intercept = intercept_sd,
                           sd_size = size_sd,
                           effect_cor = effect_cor)
      year_eff_df <- rbind(year_eff_df, tmp_df)
    } # end vital rate loop
  } # end species loop
} # end state loop

year_sds_plot <- year_eff_df %>%
  gather(effect, sdev, sd_intercept:sd_size) %>%
  mutate(state_species = paste(state,species, sep = "::"),
         vital_rate = ifelse(vital_rate=="growth","Growth","Survival")) 

ggplot(year_sds_plot, aes(x=state_species, y=sdev, color=effect))+
  geom_lollipop(point.size=2, horizontal=FALSE)+
  facet_wrap(~vital_rate)+
  scale_color_brewer(palette = "Set1", name="", labels=c("Intercept","log(size)"))+
  ylab("S.D. of Random Year Effects")+
  xlab("State::Species")+
  coord_flip()+
  theme_few()
ggsave(paste0(loc_fig_dir,"year_effects_sds.png"), width = 6, height = 3, units = "in")



####
####  CALCULATE TRANSIENT TIME TO EQUILIBRIUM AND PLOT ----
####
equil_cover <- summary_sims %>%
  dplyr::filter(simtype == "sizeXyear") %>%
  dplyr::select(state,species,mean_cover)

all_files <- list.files(paste0(results_dir,"transient_runs/"))
all_sims  <- {}
for(do_file in all_files){
  file_info <- unlist(str_split(do_file,"_"))
  state     <- file_info[1]
  species   <- file_info[2]
  simtype   <- unlist(str_split(file_info[3],"[.]"))[1]
  
  ##  Rename states with formatting for figures
  if(state == "arizona"){state <- "Arizona"}
  if(state == "kansas"){state <- "Kansas"}
  if(state == "montana"){state <- "Montana"}
  if(state == "newmexico"){state <- "New Mexico"}
  if(state == "idaho"){state <- "Idaho"}
  
  tmp <- readRDS(paste0(results_dir,"transient_runs/",do_file)) %>%
    mutate(state = state, species = species, simtype = simtype) %>%
    mutate(sim_id = rep(1:(length(iteration)/100), each = 100),
           sim_iter = rep(1:100, times = length(iteration)/100))
  
  if("cov_kids" %in% colnames(tmp)){
    tmp <- tmp %>%
      dplyr::mutate(cover = cov_kids+cov_adults,
                    abundance = N_kids+N_adults) %>%
      dplyr::select(iteration,year,cover,abundance,state,species,simtype,sim_id,sim_iter)
  }
  
  all_sims <- rbind(all_sims, tmp)
}

##  Find time at which transient cover reaches equilibrium cover
equil_time <- all_sims %>%
  filter(sim_id > 1) %>%
  dplyr::select(-abundance) %>%
  left_join(equil_cover, by = c("state", "species")) %>%
  mutate(hitit = ifelse(cover < mean_cover, "not_yet", "there")) %>%
  group_by(state,species,simtype,sim_id, occurence = hitit=="there") %>%
  arrange(state,species,simtype,sim_id) %>%
  mutate(count = row_number(), 
         first_in_group = occurence & count==1) %>%
  filter(first_in_group == TRUE) %>%
  group_by(state,species,simtype) %>%
  summarise(mean_time = mean(sim_iter),
            hi_time   = mean(sim_iter) + sd(sim_iter),
            lo_time   = mean(sim_iter) - sd(sim_iter))

transient <- ggplot(equil_time, aes(x = species, y = mean_time, fill = simtype))+
  geom_col(position = position_dodge(0.9), color = "black")+
  geom_errorbar(aes(ymin = lo_time, ymax = hi_time, group = simtype), position = position_dodge(0.9), width=0.1)+
  # scale_fill_brewer(palette = "Set1", name = "", labels=c("Size*Year", "No Size*Year"))+
  # scale_fill_manual(values=c("grey50","white"), name = "", labels=c("Size*Year", "No Size*Year"))+
  scale_fill_brewer(type = "qual", labels=c("Size*Year", "No Size*Year"), name = "")+
  ylab("Time to Equilibrium")+
  xlab("Species")+
  facet_wrap(~state, scales = "free", nrow=1)+
  theme_few()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(loc_fig_dir,"transient_times.png"), plot = transient, width = 9, height = 2.5, units = "in")


## Combine plots and save
outplot <- plot_grid(equilibrium, transient, nrow = 2, labels = "AUTO")
ggsave(paste0(ms_fig_dir,"equil_transient_combo.pdf"), plot = outplot, width = 9, height = 5, units = "in")

##  Write results to file
sink(file = paste0(ms_fig_dir,"transient_means.txt"))
timestamp()
cat("\n\n")
cat("Growth\n\n")
print(as.data.frame(equil_time))
sink()


####
####  LOOK AT SPECIES SYNCHRONY WITH AND WITHOUT SIZE BY YEAR EFFECTS ----
####
# synch_df <- {} # empty object for storage
# states <- unique(cut_sims$state)
# for(do_state in states){
#   tmp_df <- codyn::synchrony(df = filter(cut_sims, state == do_state), 
#                              time.var = "iteration", 
#                              species.var = "species", 
#                              abundance.var = "cover", 
#                              replicate.var = "simtype")
#   tmp_df$state <- do_state
#   synch_df <- rbind(synch_df, tmp_df)
# }


