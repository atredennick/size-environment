###############################################################################
##  plot_sxy_statistics.R: script to plot the results from
##  fitting growth and survival models with size*year interactions.
##  Plots focus on showing the variance for small and large plants.
##
## -----------------------------------------------------------------------------
##  Created by: Andrew Tredennick (atredenn@gmail.com)
##  Created on: December 5, 2017
################################################################################



####
####  LOAD LIBRARIES -----------------------------------------------------------
####
library(tidyverse) # Data science packages
library(dplyr)     # Data summarizing and manipulation functions
library(stringr)   # Working with strings
library(ggthemes)  # Pleasing defaults for ggplot2
library(cowplot)   # For combining ggplot2 objects
library(xtable)    # Saves dataframes as LaTeX tables



####
####  SET DIRECTORY INFORMATION ------------------------------------------------
####
main_dir    <- root
results_dir <- "../results/"
fig_dir <- "../figures/"



####
####  PLOT YEARLY ANOMALIES FOR SMALL AND LARGE PLANTS -------------------------
####
all_files <- list.files(results_dir)
do_keep   <- grep("anomalies", all_files)
all_files <- all_files[do_keep]

grow <- readRDS(paste0(results_dir,all_files[1]))[["year_only"]] %>%
  mutate(vital = "grow", model = "year_only", meanslo = 0, sdslo = 0) %>%
  dplyr::select(groupID,groupFctr,meanint,sdint,meanslo,sdslo,smallPred,
                largePred,smallPredUp,largePredUp,smallPredDn,largePredDn,
                year,state,species,vital,model) %>%
  rbind(readRDS(paste0(results_dir,all_files[1]))[["size_by_year"]] %>%
          mutate(vital = "grow", model = "size_by_year"))

surv <- readRDS(paste0(results_dir,all_files[2]))[["year_only"]] %>%
  mutate(vital = "surv", model = "year_only", meanslo = 0, sdslo = 0) %>%
  dplyr::select(groupID,groupFctr,meanint,sdint,meanslo,sdslo,smallPred,
                largePred,smallPredUp,largePredUp,smallPredDn,largePredDn,
                year,state,species,vital,model) %>%
  rbind(readRDS(paste0(results_dir,all_files[2]))[["size_by_year"]] %>%
          mutate(vital = "surv", model = "size_by_year"))

small_ones <- rbind(grow,surv) %>%
  filter(model == "size_by_year") %>%
  dplyr::select(year,state,species,vital,smallPred,smallPredUp,smallPredDn) %>%
  dplyr::rename(mean_pred = smallPred,
                upper_ci  = smallPredUp,
                lower_ci  = smallPredDn) %>%
  mutate(plant_size = "Small Plants")

large_ones <- rbind(grow,surv) %>%
  filter(model == "size_by_year") %>%
  dplyr::select(year,state,species,vital,largePred,largePredUp,largePredDn) %>%
  dplyr::rename(mean_pred = largePred,
                upper_ci  = largePredUp,
                lower_ci  = largePredDn) %>%
  mutate(plant_size = "Large Plants")

all_sims <- rbind(small_ones,large_ones) %>%
  mutate(state_spp = paste(state,species,sep="::"))

# ggplot(filter(all_sims, vital=="grow"), aes(x=year,y=mean_pred, color=plant_size))+
#   geom_hline(aes(yintercept=0))+
#   geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width=0)+
#   geom_point()+
#   facet_wrap(~state_spp, scales = "free")+
#   xlab("Year")+
#   ylab("Annual Anomaly")+
#   scale_color_manual(name=NULL,values = c("black","red"))+
#   ggtitle("Growth Anomalies")+
#   theme_few()
# ggsave(paste0(loc_fig_dir,"growth_anomalies_byyear.pdf"), width = 10, height = 8, units = "in")

# ggplot(filter(all_sims, vital=="surv"), aes(x=year,y=mean_pred, color=plant_size))+
#   geom_hline(aes(yintercept=0))+
#   geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width=0)+
#   geom_point()+
#   facet_wrap(~state_spp, scales = "free")+
#   xlab("Year")+
#   ylab("Annual Anomaly")+
#   scale_color_manual(name=NULL,values = c("black","red"))+
#   ggtitle("Survival Anomalies")+
#   theme_few()
# ggsave(paste0(loc_fig_dir,"survival_anomalies_byyear.pdf"), width = 10, height = 8, units = "in")



####
####  PLOT SD OF LARGE PLANT ANOMALIES VS. SMALL PLANTS ------------------------
####
sd_sims <- all_sims %>%
  group_by(state, species, vital, plant_size) %>%
  summarise(sd_anom = sd(mean_pred)) %>%
  spread(plant_size, sd_anom)
colnames(sd_sims) <- c("state","species","vital","sd_large","sd_small")

##  ANOVA with plant size as covariate
model_data <- sd_sims %>%
  gather(key = plant_size, value = sd_rye, sd_large:sd_small)
grow_rye_mod <- lm(sd_rye ~ plant_size, data = filter(model_data,vital=="grow"))
surv_rye_mod <- lm(sd_rye ~ plant_size, data = filter(model_data,vital=="surv"))

##  Write ANOVA results to file
sink(file = paste0(results_dir,"rye_anovas.txt"))
timestamp()
cat("\n\n")
cat("Growth\n\n")
print(anova(grow_rye_mod))
cat("\n\n\n")
cat("Survival\n\n")
print(anova(surv_rye_mod))
sink()


##  Plots
grow_sds <- ggplot(data = filter(sd_sims,vital=="grow"), aes(x=sd_large, y=sd_small, fill=state))+
  geom_abline(aes(intercept=0, slope=1), lty=2)+
  geom_point(size=3, shape=21)+
  annotate("text", x = Inf, y = -Inf, label = "S.D. ~ size: P = 0.627",size = 3,hjust=1.1, vjust=-1)+
  scale_y_continuous(limits=c(0,1))+
  scale_x_continuous(limits=c(0,1))+
  # scale_fill_brewer(palette = "Set1")+
  scale_fill_colorblind(name=NULL)+
  xlab("S.D. of Large Plant Anomalies")+
  ylab("S.D. of Small Plant Anomalies")+
  ggtitle("Growth")+
  theme_few()+
  guides(fill=guide_legend(override.aes = list(size=2), keywidth=0.12,
                           keyheight=0.12,
                           default.unit="inch"))+
  theme(legend.position = c(0.15,0.86),
        legend.text=element_text(size=8))

surv_sds <- ggplot(filter(sd_sims,vital=="surv"), aes(x=sd_large, y=sd_small, fill=state))+
  geom_abline(aes(intercept=0, slope=1), lty=2)+
  geom_point(size=3, shape=21)+
  annotate("text", x = Inf, y = -Inf, label = "S.D. ~ size: P = 0.002",size = 3,hjust=1.1, vjust=-1)+
  scale_y_continuous(limits=c(0,0.25))+
  scale_x_continuous(limits=c(0,0.25))+
  # scale_fill_brewer(palette = "Set1")+
  scale_fill_colorblind()+
  xlab("S.D. of Large Plant Anomalies")+
  ylab("S.D. of Small Plant Anomalies")+
  guides(fill=FALSE)+
  ggtitle("Survival")+
  theme_few()

plot_grid(grow_sds,surv_sds,labels = "AUTO")
ggsave(paste0(fig_dir,"sd_anomalies_all.pdf"), width = 8, height = 4, units = "in")



####
####  SCATTERPLOT OF ANOMALIES OF SMALL AND LARGE PLANTS FOR EACH SPECIES ------
####
small_large_compare <- all_sims %>%
  dplyr::select(-upper_ci, -lower_ci) %>%
  spread(plant_size, mean_pred)

vital_compare <- all_sims %>%
  dplyr::select(-upper_ci, -lower_ci) %>%
  spread(vital, mean_pred)

mycols <- c("#ff8c77",
            "#271f8b",
            "#96c63f",
            "#b45b7f",
            "#fece59",
            "#47245f",
            "#007519",
            "#eb408a",
            "#61ee90",
            "#601900",
            "#7c93ff",
            "#e1d681",
            "#80188d",
            "#00bf8d",
            "#ce46ae",
            "#7a9200",
            "#d16ee1",
            "#fe7846",
            "#005aa5")
growth_plots <- ggplot(filter(small_large_compare, vital=="grow"), aes(x=`Small Plants`, y=`Large Plants`, fill=species))+
  geom_hline(aes(yintercept=0), linetype="dashed", lwd=0.25)+
  geom_vline(aes(xintercept=0), linetype="dashed", lwd=0.25)+
  geom_point(color="grey25", shape=21, size=2)+
  scale_color_manual(values = mycols)+
  facet_wrap(~state, scales = "free", nrow=1)+
  guides(fill=FALSE)+
  ggtitle("Growth")+
  theme_few()+
  theme(axis.text=element_text(size=7))

surv_plots <- ggplot(filter(small_large_compare, vital=="surv"), aes(x=`Small Plants`, y=`Large Plants`, fill=species))+
  geom_hline(aes(yintercept=0), linetype="dashed", lwd=0.25)+
  geom_vline(aes(xintercept=0), linetype="dashed", lwd=0.25)+
  geom_point(color="grey25", shape=21, size=2)+
  scale_color_manual(values = mycols)+
  facet_wrap(~state, scales = "free", nrow=1)+
  guides(fill=FALSE)+
  ggtitle("Survival")+
  theme_few()+
  theme(axis.text=element_text(size=7))

large_plots <- ggplot(filter(vital_compare, plant_size=="Large Plants"), aes(x=grow, y=surv, fill=species))+
  geom_hline(aes(yintercept=0), linetype="dashed", lwd=0.25)+
  geom_vline(aes(xintercept=0), linetype="dashed", lwd=0.25)+
  geom_point(color="grey25", shape=21, size=2)+
  scale_color_manual(values = mycols)+
  facet_wrap(~state, scales = "free", nrow=1)+
  guides(fill=FALSE)+
  xlab("Growth")+
  ylab("Survival")+
  ggtitle("Large Plants")+
  theme_few()+
  theme(axis.text=element_text(size=7))

small_plots <- ggplot(filter(vital_compare, plant_size=="Small Plants"), aes(x=grow, y=surv, fill=species))+
  geom_hline(aes(yintercept=0), linetype="dashed", lwd=0.25)+
  geom_vline(aes(xintercept=0), linetype="dashed", lwd=0.25)+
  geom_point(color="grey25", shape=21, size=2)+
  scale_color_manual(values = mycols)+
  facet_wrap(~state, scales = "free", nrow=1)+
  guides(fill=FALSE)+
  xlab("Growth")+
  ylab("Survival")+
  ggtitle("Small Plants")+
  theme_few()+
  theme(axis.text=element_text(size=7))

suppressWarnings( # ignore warning about a few non-aligning years for growth and survival
  # out_plot <- plot_grid(growth_plots,surv_plots,small_plots,large_plots,nrow = 4,labels = "AUTO")
  out_plot <- plot_grid(growth_plots,surv_plots,nrow = 2,labels = "AUTO")
)
ggsave(paste0(fig_dir,"small_large_corrs.pdf"),plot = out_plot, width = 10, height = 5, units = "in")

# test_data <- filter(small_large_compare, state=="NewMexico", vital=="grow")
# t.test(test_data$`Large Plants`, test_data$`Small Plants`, paired=TRUE)



####
####  CALCULATE CORRELATIONS OF RANDOM YEAR ANAMOLIES --------------------------
####
plant_size_corrs <- small_large_compare %>%
  group_by(state,species,vital) %>%
  summarise(anom_cor = cor(`Large Plants`,`Small Plants`, 
                           use = "pairwise.complete.obs", 
                           method = "pearson"),
            cor_pvalue = cor.test(`Large Plants`,`Small Plants`, 
                                  use = "pairwise.complete.obs", 
                                  method = "pearson")$p.value)

vital_rate_corrs <- vital_compare %>%
  group_by(state,species,plant_size) %>%
  summarise(anom_cor = cor(grow,surv,
                           use = "pairwise.complete.obs",
                           method = "pearson"),
            cor_pvalue = cor.test(grow,surv,
                                  use = "pairwise.complete.obs",
                                  method = "pearson")$p.value)

##  Write to text file
sink(file = paste0(results_dir,"rye_corrs.txt"))
timestamp()
cat("\n\n")
cat("Among sizes\n\n")
print(as.data.frame(plant_size_corrs))
cat("\n\n\n")
cat("Among vital rates\n\n")
print(as.data.frame(vital_rate_corrs))
sink()


####
####  TALLY NUMBER OF CONSISTENT AND INCONSISTENT RANDOM YEAR ANAMOLIES --------
####
small_large_signs <- small_large_compare %>%
  dplyr::mutate(large_sign = sign(`Large Plants`),
                small_sign = sign(`Small Plants`)) %>%
  dplyr::mutate(sign_add = large_sign+small_sign) %>% 
  group_by(state,species,vital) %>%
  summarise(offs = length(which(sign_add == 0)),
            ons = length(which(sign_add != 0))) %>%
  dplyr::mutate(off_perc = offs / (offs+ons) * 100,
                on_perc = ons / (offs+ons) * 100)

signs_by_state <- small_large_compare %>%
  dplyr::mutate(large_sign = sign(`Large Plants`),
                small_sign = sign(`Small Plants`)) %>%
  dplyr::mutate(sign_add = large_sign+small_sign) %>% 
  group_by(state,vital) %>%
  summarise(offs = length(which(sign_add == 0)),
            ons = length(which(sign_add != 0))) %>%
  dplyr::mutate(off_perc = offs / (offs+ons) * 100,
                on_perc = ons / (offs+ons) * 100)

signs_by_vital <- small_large_compare %>%
  dplyr::mutate(large_sign = sign(`Large Plants`),
                small_sign = sign(`Small Plants`)) %>%
  dplyr::mutate(sign_add = large_sign+small_sign) %>% 
  group_by(vital) %>%
  summarise(offs = length(which(sign_add == 0)),
            ons = length(which(sign_add != 0))) %>%
  dplyr::mutate(off_perc = offs / (offs+ons) * 100,
                on_perc = ons / (offs+ons) * 100)

sink(file = paste0(results_dir,"rye_diag_tallies.txt"))
timestamp()
cat("\n\n")
cat("Grouped by state, species, and vital rate\n\n")
print(as.data.frame(small_large_signs))
cat("\n\n\n")
cat("Grouped by state and vital rate\n\n")
print(signs_by_state)
cat("\n\n\n")
cat("Grouped by vital rate\n\n")
print(signs_by_vital)
sink()


surv_grow_signs <- vital_compare %>%
  dplyr::mutate(grow_sign = sign(`grow`),
                surv_sign = sign(`surv`)) %>%
  dplyr::mutate(sign_add = grow_sign+surv_sign) %>% 
  group_by(state,species,plant_size) %>%
  summarise(offs = length(which(sign_add == 0)),
            ons = length(which(sign_add != 0))) %>%
  dplyr::mutate(off_perc = offs / (offs+ons) * 100,
                on_perc = ons / (offs+ons) * 100)

surv_grow_signs_by_state <- vital_compare %>%
  dplyr::mutate(grow_sign = sign(`grow`),
                surv_sign = sign(`surv`)) %>%
  dplyr::mutate(sign_add = grow_sign+surv_sign) %>% 
  group_by(state,plant_size) %>%
  summarise(offs = length(which(sign_add == 0)),
            ons = length(which(sign_add != 0))) %>%
  dplyr::mutate(off_perc = offs / (offs+ons) * 100,
                on_perc = ons / (offs+ons) * 100)

surv_grow_signs_only <- vital_compare %>%
  dplyr::mutate(grow_sign = sign(`grow`),
                surv_sign = sign(`surv`)) %>%
  dplyr::mutate(sign_add = grow_sign+surv_sign) %>% 
  group_by(plant_size) %>%
  summarise(offs = length(which(sign_add == 0)),
            ons = length(which(sign_add != 0))) %>%
  dplyr::mutate(off_perc = offs / (offs+ons) * 100,
                on_perc = ons / (offs+ons) * 100)

sink(file = paste0(results_dir,"rye_diag_tallies_within_sizes.txt"))
timestamp()
cat("\n\n")
cat("Grouped by state, species, and vital rate\n\n")
print(as.data.frame(surv_grow_signs))
cat("\n\n\n")
cat("Grouped by state and vital rate\n\n")
print(surv_grow_signs_by_state)
cat("\n\n\n")
cat("Grouped by vital rate\n\n")
print(surv_grow_signs_only)
sink()



### SI PLOTS ###

####
####  COMPARE ANOMALIES FROM BOTH MODELS ---------------------------------------
####

##  Large plants
large_anoms <- rbind(grow,surv) %>%
  dplyr::select(year,state,species,vital,model,largePred,largePredUp,largePredDn) %>%
  dplyr::rename(mean_pred = largePred,
                upper_ci  = largePredUp,
                lower_ci  = largePredDn) %>%
  mutate(plant_size = "Large plants") %>%
  group_by(state,species,vital,model,plant_size) %>%
  summarise(sd_anom = sd(mean_pred)) %>%
  spread(model, sd_anom) %>%
  ungroup() %>%
  mutate(vital = ifelse(vital == "grow", "Growth", "Survival"))

##  Small plants
small_anoms <- rbind(grow,surv) %>%
  dplyr::select(year,state,species,vital,model,smallPred,smallPredUp,smallPredDn) %>%
  dplyr::rename(mean_pred = smallPred,
                upper_ci  = smallPredUp,
                lower_ci  = smallPredDn) %>%
  mutate(plant_size = "Small plants") %>%
  group_by(state,species,vital,model,plant_size) %>%
  summarise(sd_anom = sd(mean_pred)) %>%
  spread(model, sd_anom) %>%
  ungroup() %>%
  mutate(vital = ifelse(vital == "grow", "Growth", "Survival"))

all_sdanoms <- rbind(large_anoms, small_anoms) %>%
  dplyr::mutate(state = ifelse(state == "NewMexico", "New Mexico", as.character(state)))

mod_comp_anoms <- ggplot(all_sdanoms, aes(x = year_only, y = size_by_year, fill=state))+
  geom_abline(aes(intercept = 0, slope = 1), linetype="dashed", lwd=0.25)+
  geom_point(color="grey25", shape=21, size=2)+
  facet_grid(plant_size~vital, scales = "free")+
  xlab("S.D. of anomalies from year only model")+
  ylab("S.D. of anomalies from size*year model")+
  scale_fill_brewer(palette = "Set1", name = NULL)+
  theme_few()+
  theme(legend.position = c(0.1,0.9),
        legend.background = element_rect(colour = NA, fill = NA),
        legend.text=element_text(size=7),
        legend.key.size = unit(0.4, "cm"))
# ggsave(paste0(ms_fig_dir,"model_comp_anoms.pdf"), width = 6, height = 5, units = "in")


##  Plot the SDs for size-independent models
all_sdanoms <- all_sdanoms %>%
  mutate(state_spp = paste(state, species, sep = ", "))
ggplot(all_sdanoms, aes(x = state_spp, y = year_only, fill = state))+
  geom_col()+
  facet_grid(plant_size~vital)+
  scale_fill_brewer(palette = "Set1", name = NULL)+
  coord_flip()+
  guides(fill = FALSE)+
  xlab(NULL)+
  ylab("S.D. of Anomaly")+
  theme_few()
ggsave(paste0(ms_fig_dir,"sizeIndependent_sd_anomalies.pdf"), height = 6, width = 6, units = "in")

####
####  CALCULATE CORRELATION OF THE ANOMALIES BETWEEN MODELS --------------------
####
large_anom_cors <- rbind(grow, surv) %>%
  dplyr::select(year,state,species,vital,model,largePred) %>%
  spread(model, largePred) %>%
  group_by(state, species, vital) %>%
  summarise(anom_cor = cor(size_by_year, year_only)) %>%
  ungroup() %>%
  dplyr::mutate(plant_size = "Large plants")

small_anom_cors <- rbind(grow, surv) %>%
  dplyr::select(year,state,species,vital,model,smallPred) %>%
  spread(model, smallPred) %>%
  group_by(state, species, vital) %>%
  summarise(anom_cor = cor(size_by_year, year_only)) %>%
  ungroup() %>%
  dplyr::mutate(plant_size = "Small plants")

all_anom_cors <- rbind(large_anom_cors, small_anom_cors) %>%
  dplyr::mutate(vital = ifelse(vital == "grow", "Growth", "Survival"),
                state = ifelse(state == "NewMexico", "New Mexico", as.character(state)),
                state_spp = paste(state, species, sep = ", "))

anom_corplot <- ggplot(all_anom_cors, aes(x = state_spp, y = anom_cor, fill = state))+
  geom_col()+
  geom_hline(aes(yintercept = 0))+
  geom_hline(aes(yintercept = 0.5), linetype = "dashed", color = "grey25")+
  xlab("")+
  ylab("Correlation of anomalies")+
  scale_fill_brewer(palette = "Set1", name = NULL)+
  scale_y_continuous(breaks = c(0,0.5,1))+
  facet_grid(plant_size~vital)+
  coord_flip()+
  guides(fill = FALSE)+
  theme_few()
# ggsave(paste0(ms_fig_dir,"model_cor_anoms.pdf"), plot = anom_corplot, width = 5, height = 6, units = "in")

##  Combine model-type comparison anomaly plots
theplots <- list(anom_corplot, NULL, mod_comp_anoms)
mod_anom_panel <- plot_grid(plotlist = theplots, 
                            ncol = 3,
                            labels = c("A","","B"),
                            rel_widths = c(1,0.1,1))
ggsave(paste0(ms_fig_dir,"model_cor_anoms.pdf"), 
       plot = mod_anom_panel, 
       width = 11, 
       height = 5, 
       units = "in")



