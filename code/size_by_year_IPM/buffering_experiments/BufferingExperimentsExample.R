################################################################################
##  BufferingExperimentsExample.R: R script to produce example plots of the
##  buffering experiments we performed.
##
## -----------------------------------------------------------------------------
##  Author: Stephen Ellner
##          Andrew Tredennick (atredenn@gmail.com)
##  Date created: December 1, 2017
################################################################################



####
####  LOAD LIBRARIES -----------------------------------------------------------
####
library(tidyverse)
library(dplyr)
library(cowplot)
library(ggthemes)



####
####  SET DIRECTORIES ----------------------------------------------------------
####
root    <- ifelse(.Platform$OS.type=="windows","c:/repos","~/repos")
fig_dir <- paste0(root,"/drivers/Manuscripts/SizebyClimate/figures/")
setwd(fig_dir) # set working directory



####
####  DEFINE SOME VARIABLES ----------------------------------------------------
####
##  Plant sizes 
z = runif(1000,-2,5)
z = sort(z)

##  Baseline coefficients 
a0 = 1
a1 = 0.5

##  Random year-effects on intercept 
Nyears = 50
alpha = qnorm(seq(.05,.95,length=Nyears))

##  Baseline mean square (MS) anomaly 
Ab = sqrt(mean(alpha^2))

##  Parameters for defining SxY models 
zL = quantile(z,0.05) # small plant size
zU = quantile(z,0.95) # large plant size
c = 3/(zU-zL)
zbar = ztilde = (zL + zU)/2



####
####  CALCULATE BASELINE MODEL PREDICTIONS -------------------------------------
####
yhat = matrix(NA,1000,Nyears) 
for(j in 1:Nyears) yhat[,j] = a0 + a1*z + alpha[j]

saved_outs <- as.data.frame(yhat)
colnames(saved_outs) <- as.character(1:ncol(saved_outs))
saved_outs <- saved_outs %>%
  dplyr::mutate(z_t0 = z,
                `13` = a0+a1*z) %>%
  gather(key = year, value = yhat, -z_t0) %>%
  dplyr::mutate(year_group = ifelse(year == "13", "mean", "random"),
                rand_eff = "Baseline")


####
####  SMALL PLANTS MORE VARIABLE -----------------------------------------------
####
##  Compute anomalies for q=1
anomalies = numeric(0)
for(j in 1:Nyears) {
    aj = alpha[j]*exp(-c*(z-ztilde)) 
    anomalies = c(anomalies,aj)
}   
A = sqrt(mean(anomalies^2))
q = Ab/A; # scale up to match baseline RMS anomaly. 

yhat = matrix(NA,1000,Nyears)
for(j in 1:Nyears) yhat[,j]=a0 + a1*z + q*alpha[j]*exp(-c*(z-ztilde)) 

tmp_outs <- as.data.frame(yhat)
colnames(tmp_outs) <- as.character(1:ncol(tmp_outs))
tmp_outs <- tmp_outs %>%
  dplyr::mutate(z_t0 = z,
                `13` = a0+a1*z) %>%
  gather(key = year, value = yhat, -z_t0) %>%
  dplyr::mutate(year_group = ifelse(year == "13", "mean", "random"),
                rand_eff = "Small plants variable")

saved_outs <- rbind(saved_outs, tmp_outs)



####
####  LARGE PLANTS MORE VARIABLE -----------------------------------------------
####
##  Compute anomalies for q=1
anomalies = numeric(0); 
for(j in 1:Nyears) {
    aj = alpha[j]*exp(c*(z-ztilde)) 
    anomalies=c(anomalies,aj)
}   
A = sqrt(mean(anomalies^2)) 
q = Ab/A # scale up to match baseline RMS anomaly. 

yhat=matrix(NA,1000,Nyears) 
for(j in 1:Nyears) yhat[,j]=a0 + a1*z + q*alpha[j]*exp(c*(z-ztilde)) 

tmp_outs <- as.data.frame(yhat)
colnames(tmp_outs) <- as.character(1:ncol(tmp_outs))
tmp_outs <- tmp_outs %>%
  dplyr::mutate(z_t0 = z,
                `13` = a0+a1*z) %>%
  gather(key = year, value = yhat, -z_t0) %>%
  dplyr::mutate(year_group = ifelse(year == "13", "mean", "random"),
                rand_eff = "Large plants variable")

saved_outs <- rbind(saved_outs, tmp_outs)



####
####  SIGN FLIP ----------------------------------------------------------------
####
zbar = (zL+zU)/2
dz = (zU-zL)

##  Compute anomalies for q=1
anomalies=numeric(0)
for(j in 1:Nyears) {
    aj = 2*alpha[j]*(z-zbar)/dz 
    anomalies = c(anomalies,aj)
}  
A = sqrt(mean(anomalies^2))
q = Ab/A # scale up to match baseline RMS anomaly. 

yhat=matrix(NA,1000,Nyears)
for(j in 1:Nyears) yhat[,j]=a0 + a1*z + q*2*alpha[j]*(z-zbar)/dz 

tmp_outs <- as.data.frame(yhat)
colnames(tmp_outs) <- as.character(1:ncol(tmp_outs))
tmp_outs <- tmp_outs %>%
  dplyr::mutate(z_t0 = z,
                `13` = a0+a1*z) %>%
  gather(key = year, value = yhat, -z_t0) %>%
  dplyr::mutate(year_group = ifelse(year == "13", "mean", "random"),
                rand_eff = "Sign flip")

saved_outs <- rbind(saved_outs, tmp_outs)


saved_outs <- saved_outs %>%
  filter(rand_eff != "Sign flip")


####
####  PLOT THE RESPONSES ON LOG AND ARITHMETIC SCALES --------------------------
####
logged <- ggplot(saved_outs, aes(x = z_t0, y = yhat, group = year, color = year_group))+
  geom_line(aes(linetype = year_group, size = year_group))+
  facet_wrap(~rand_eff, nrow = 1, scales = "free")+
  guides(color = FALSE, linetype = FALSE, size = FALSE)+
  scale_size_manual(values = c(1,0.5))+
  scale_linetype_manual(values = c(1,3))+
  xlab(expression(paste("Size, ",log(z[t]))))+
  ylab(expression(paste("Size, ",log(z[t+1]))))+
  theme_few()

arith <- ggplot(saved_outs, aes(x = exp(z_t0), y = exp(yhat), group = year, color = year_group))+
  geom_line(aes(linetype = year_group, size = year_group))+
  facet_wrap(~rand_eff, nrow = 1, scales = "free")+
  guides(color = FALSE, linetype = FALSE, size = FALSE)+
  scale_size_manual(values = c(1,0.5))+
  scale_linetype_manual(values = c(1,3))+
  xlab(expression(paste("Size, ",z[t])))+
  ylab(expression(paste("Size, ",z[t+1])))+
  theme_few()

# outplot <- plot_grid(logged, arith, nrow = 2, labels = "AUTO", align = "v")
# ggsave("./buffering_examples.pdf", plot = outplot, width = 8.5, height = 5)



####
####  LOOK AT PLOTS OF LOGIT TRANSFORM -----------------------------------------
####
antilogit <- function(x) { exp(x) / (1 + exp(x) ) }
saved_outs <- saved_outs %>%
  mutate(logit_yhat = antilogit(yhat))

surv <- ggplot(saved_outs, aes(x = z_t0, y = logit_yhat, group = year, color = year_group))+
  geom_line(aes(linetype = year_group, size = year_group))+
  facet_wrap(~rand_eff, nrow = 1, scales = "free")+
  guides(color = FALSE, linetype = FALSE, size = FALSE)+
  scale_size_manual(values = c(1,0.5))+
  scale_linetype_manual(values = c(1,3))+
  xlab(expression(paste("Size, ",log(z[t]))))+
  ylab("Probability of survival")+
  theme_few()

# outplot <- plot_grid(logged, arith, surv, nrow = 3, labels = "AUTO", align = "v")
# ggsave("./buffering_examples.pdf", plot = outplot, width = 8.5, height = 7.5)



####
####  PLOT THE DISTRIBUTION OF Z_HAT AT SMALL AND LARGE SLICE ------------------
####
# Get the closest real values to the quantile estimates for large and small
zL_real <- saved_outs$z_t0[which.min(abs(saved_outs$z_t0 - zL))]
zU_real <- saved_outs$z_t0[which.min(abs(saved_outs$z_t0 - zU))] 

# Subset out only large and small plants
small_large_sub <- saved_outs %>%
  filter(z_t0 %in% c(zL_real, zU_real)) %>%
  filter(year_group == "random") %>%
  mutate(size = ifelse(z_t0 == zL_real, "Small plant,", "Large plant,"))

# Plot their yhat distributions
ggplot(small_large_sub, aes(x = exp(yhat), fill = size))+
  geom_histogram(color = "white", bins = 30)+
  facet_wrap(size~rand_eff, scales = "free")+
  scale_fill_brewer(type = "qual")+
  theme_few()+
  guides(fill = FALSE)+
  xlab(expression(paste("Predicted size, ",z[t+1])))+
  ylab("Frequency")



