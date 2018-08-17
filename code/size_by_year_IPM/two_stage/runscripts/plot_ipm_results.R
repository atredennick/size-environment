##  plot_ipm_results.R: script to plot results from IPM simulations with
##  and without size-by-year effects.



##  Clear everything
rm(list = ls(all.names = TRUE))


####
####  LOAD LIBRARIES ----
####
library(tidyverse)
library(dplyr)
library(stringr)
library(ggthemes)
library(ggjoy)



####
####  SET DIRECTORIES ----
####
root        <- ifelse(.Platform$OS.type=="windows","c:/repos","~/repos")
work_dir    <- paste0(root,"/drivers/empirical/size_by_year_IPM/two_stage/")
results_dir <- paste0(root,"/drivers/empirical/size_by_year_IPM/two_stage/ipm_results/")
out_dir     <- paste0(root,"/drivers/empirical/size_by_year_IPM/two_stage/figures/")

setwd(work_dir) # set working directory
dir.create(file.path(out_dir), showWarnings = FALSE) # create a folder for figures



####
####  READ IN AND COMBINE RESULTS FOR EACH SPECIES AND STATE ----
####
all_files <- list.files(results_dir)
all_sims  <- {}
for(do_file in all_files){
  file_info <- unlist(str_split(do_file,"_"))
  state     <- file_info[1]
  species   <- file_info[2]
  simtype   <- unlist(str_split(file_info[3],"[.]"))[1]
  
  tmp <- readRDS(paste0(results_dir,do_file)) %>%
    mutate(state = state, species = species, simtype = simtype)
  
  all_sims <- rbind(all_sims, tmp)
}
cut_sims <- filter(all_sims, iteration>199) %>%
  dplyr::mutate(cover = cov_kids+cov_adults)



####
####  CALCULATE STD. DEV. OF EACH RUN ----
####
summary_sims <- cut_sims %>%
  group_by(state, species, simtype) %>%
  summarise(sd_cover = sd(cover)) %>%
  mutate(sim_text = ifelse(simtype=="sizeXyear", "with SxY", "without SxY")) %>%
  mutate(sd_text = paste("S.D.", sim_text, "=",round(sd_cover,3))) %>%
  ungroup()
summary_sims$xpos <- c(8, 8, 13, 13, 5.5, 5.5, 5,5)
summary_sims$ypos <- c(600,550,1000,920,800,730,600,550)


####
####  MAKE THE PLOTS ----
####
mycols <- c("#ff000050", "#0000ff50")

ggplot(cut_sims, aes(x=cover*100, fill=simtype))+
  geom_histogram(alpha=0.4, color="black", size=0.1, bins=20)+
  geom_text(data = summary_sims, aes(x=xpos, y=ypos, label = sd_text), size = 2.5)+
  facet_wrap(~species, scales = "free")+
  scale_fill_manual(values = mycols, name=NULL, labels=c("With SxY", "Without SxY"))+
  xlab("Cover (%)")+
  ylab("Frequency")+
  theme_few()+
  theme(legend.position = "top",
        legend.key.size = unit(0.2, "cm"))
#ggsave(paste0(out_dir,"idaho_one-stage_histograms.pdf"), width = 6, height = 5, units = "in")


ggplot(cut_sims, aes(x=species,y=cover*100,fill=simtype))+
  geom_boxplot(width=0.5)+
  scale_fill_brewer(palette = "Set1")+
  ylab("Cover (%)")+
  xlab("Species")+
  ggtitle("Idaho")
ggsave(paste0(out_dir,"idaho_two-stage_boxplots.pdf"), width = 6, height = 5, units = "in")

##  Example time series from IPM
ggplot(filter(cut_sims, species=="ARTR"), aes(x=iteration, y=cover, color=simtype))+
  geom_line()



####
####  LOOK AT VARIANCE OF RANDOM YEAR EFFECTS -----
####
all_param_files <- as.data.frame(list.files("./vital_rate_params/")) %>%
  dplyr::rename(fname = !!names(.[1])) %>%
  separate(fname, c("vital_rate","app","state","spp_ext"), "_") %>%
  separate(spp_ext, c("species","ext"), "[.]")

states <- unique(all_param_files$state)
year_eff_df <- {} # empty object
set_graph_pars <- function(ptype = "panel4") {
  mgp <- c(2.5, 1, 0); mar <- c(4, 4, 1.5, 1); oma <- c(0, 0, 0, 0)
  switch(ptype,
         panel4 = par(mfrow=c(2,2), mgp=mgp, mar=mar, pty="s",
                      oma=oma, bty="L", cex.lab =1.2),
         panel2 = par(mfrow=c(1,2), mgp=mgp, mar=mar, pty="s",
                      oma=oma, bty="L", cex.axis=0.85),
         panel1 = par(mfrow=c(1,1), mgp=mgp, mar=mar, pty="s",
                      oma=oma, bty="L", cex.axis=0.85))
}

pdf(file = paste0(out_dir,"idaho_random_effects_corrs.pdf"), height = 5)
set_graph_pars()
for(do_state in states){
  species <- unique(filter(all_param_files, state==do_state)$species)
  
  for(do_species in species){
    vitals <- unique(filter(all_param_files, state==do_state & species==do_species))$vital_rate
    
    for(do_vital in "growth"){
      file_parts <- filter(all_param_files, state==do_state & species==do_species & vital_rate==do_vital)
      file_name <- file_parts %>%
        unite(fname1, vital_rate:species, sep="_") %>%
        unite(fname, fname1:ext, sep=".")
      tmp_params <- readRDS(paste0("./vital_rate_params/",file_name))[["sizeXyear"]]
      intercept_sd <- sd(tmp_params$intercept)
      size_sd <- sd(tmp_params$logarea.t0)
      effect_cor <- cor(tmp_params$intercept, tmp_params$logarea.t0)
      plot(tmp_params$intercept, tmp_params$logarea.t0, xlab="Intercept", ylab="logarea.t0",main=do_species)
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
dev.off()

year_eff_df %>%
  gather(effect, sd, sd_intercept:sd_size) %>%
  ggplot(aes(x=species, y=sd, fill=effect))+
    geom_col(position = position_dodge(0.9))+
    scale_fill_brewer(palette = "Set1", name="", labels=c("Intercept","Size"))+
    xlab("Species")+
    ylab("Standard Deviation Across Years")+
    ggtitle("Idaho")+
    coord_flip()
ggsave(paste0(out_dir,"idaho_random_effects_df.pdf"), width = 6, height = 3, units = "in")


