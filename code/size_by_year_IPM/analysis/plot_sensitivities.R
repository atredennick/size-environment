################################################################################
##  plot_sensitivities.R: script to plot results from IPM simulations with
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
library(cowplot)   # Combining ggplot2 objects



####
####  SET DIRECTORIES ----
####
root        <- ifelse(.Platform$OS.type=="windows","c:/repos","~/repos")
work_dir    <- paste0(root,"/drivers/empirical/size_by_year_IPM/analysis/")
results_dir <- paste0(root,"/drivers/empirical/size_by_year_IPM/ipm_results/sensitivities/")
loc_fig_dir <- "/Users/atredenn/Google_Drive/size_year_paper/figures/" # update as needed, storing off repo
ms_fig_dir  <- paste0(root,"/drivers/Manuscripts/SizebyClimate/figures/")


setwd(work_dir) # set working directory
dir.create(file.path(loc_fig_dir), showWarnings = FALSE) # create a folder for figures
# dir.create(file.path(ms_fig_dir), showWarnings = FALSE) # create a folder for figures



####
####  DEFINE PLOTTING FUNCTION ----
####
matrix.image <- function(A, x=NULL, y=NULL, col=rainbow(100,start=0.67,end=0),
                           bw=FALSE, do.contour=FALSE, do.legend=TRUE,...) {
  if(do.legend) layout(mat=cbind(matrix(1,5,5),rep(2,5)));
  par(mar=c(6,5,3,2));
  if(is.null(x)) x=1:ncol(A);
  if(is.null(y)) y=1:nrow(A);
  nx=length(x); ny=length(y);
  x1=c(1.5*x[1]-0.5*x[2],1.5*x[nx]-0.5*x[nx-1]);
  y1=c(1.5*y[1]-0.5*y[2],1.5*y[ny]-0.5*y[ny-1]);
  if(bw) col=grey( (200:50)/200 );
  image(list(x=x,y=y,z=A),xlim=x1,ylim=rev(y1),col=col,cex.axis=1,cex.lab=1,bty="u",...);
  abline(v=range(x1)); abline(h=range(y1));
  if(do.contour) {
    labs <- round(tapply(as.numeric(A), cut(as.numeric(A), 5), max),1)
    contour(x,y,A,nlevels=6,labcex=0.5,add=TRUE);
  }
  
  if(do.legend) {
    l.y=seq(min(A),max(A),length=100);
    par(mar=c(6,2,3,1))
    image(list(x=1:2,y=l.y,z=rbind(l.y,l.y)),col=col,bty="o",xaxt="n",yaxt="n");
    axis(side=2,cex.axis=1,at=pretty(seq(min(A),max(A),length=10)));
  }
}


####
####  LOOP THROUGH FILES AND PLOT SENSITIVITIES ----
####
# pdf(file = "./figures/all_sensitivities.pdf", width = 11, height = 11)
# par(mfrow=c(4,4),tcl=-0.2,mgp=c(2,0.5,0),mar=c(3,4,1,1),cex.lab=1.2, las=1)
# all_files <- list.files(paste0(results_dir))
# for(do_file in all_files){
#   file_info <- unlist(str_split(do_file,"_"))
#   state     <- file_info[1]
#   species   <- file_info[2]
#   
#   ##  Rename states with formatting for figures
#   if(state == "arizona")   {state <- "Arizona"}
#   if(state == "kansas")    {state <- "Kansas"}
#   if(state == "montana")   {state <- "Montana"}
#   if(state == "newmexico") {state <- "New Mexico"}
#   if(state == "idaho")     {state <- "Idaho"}
#   
#   tmp_sens <- readRDS(paste0(results_dir,do_file))[["sensitivity"]]
#   tmp_u    <- readRDS(paste0(results_dir,do_file))[["size_bins"]]
#   
#   if(dim(tmp_sens)[1] != length(tmp_u)) tmp_sens <- sqrt(tmp_sens[2:101,2:101])
#   matrix.image(t(tmp_sens), tmp_u, tmp_u, do.contour = TRUE, bw=TRUE, 
#                  xlab = "Size (t), z", ylab = "Size (t+1), z\'", do.legend = FALSE)
#   title(paste(state, species, sep=", "))
# }
# dev.off()



####
####  LOOP THROUGH FILES AND PLOT ELASTICITIES ----
####
outfile <- paste0(loc_fig_dir, "all_elasticities.pdf")
pdf(file = outfile, width = 4, height = 4, onefile = TRUE)
par(mfrow=c(1,1),tcl=-0.2,mgp=c(2,0.5,0),mar=c(0,0,0,0),cex.lab=1.2, las=1)
all_files <- list.files(paste0(results_dir))
for(do_file in all_files){
  file_info <- unlist(str_split(do_file,"_"))
  state     <- file_info[1]
  species   <- file_info[2]
  
  ##  Rename states with formatting for figures
  if(state == "arizona")   {state <- "Arizona"}
  if(state == "kansas")    {state <- "Kansas"}
  if(state == "montana")   {state <- "Montana"}
  if(state == "newmexico") {state <- "New Mexico"}
  if(state == "idaho")     {state <- "Idaho"}
  
  tmp_elas <- readRDS(paste0(results_dir,do_file))[["elasticity"]]
  tmp_u    <- readRDS(paste0(results_dir,do_file))[["size_bins"]]
  
  if(dim(tmp_elas)[1] != length(tmp_u)) tmp_elas <- tmp_elas[2:101,2:101]
  
  # elas_df <- as.data.frame(t(tmp_elas))
  # colnames(elas_df) <- as.character(tmp_u)
  # elas_df <- elas_df %>%
  #   dplyr::mutate(size_t0 = tmp_u) %>%
  #   gather(key = size_t1, value = elasticity, -size_t0) %>%
  #   dplyr::mutate(size_t1 = as.numeric(size_t1))
  # ggplot(elas_df, aes(x = size_t0, y = size_t1))+
  #   geom_tile(aes(fill = elasticity))+
  #   # geom_contour(aes(z = elasticity), bins = 5)+
  #   scale_y_reverse()+
  #   scale_fill_continuous(trans = "sqrt")
  #   guides(fill = FALSE)
  matrix.image(t(sqrt(tmp_elas)), tmp_u, tmp_u, do.contour = FALSE, bw=FALSE, 
               xlab = "Size (t), z", ylab = "Size (t+1), z\'", do.legend = FALSE, las =1)
  title(paste(state, species, sep=", "))
  contour(tmp_u,tmp_u,t(tmp_elas),nlevels=6,labcex=1,add=TRUE)    
}
dev.off()


# 
# ####
# ####  LOOP THROUGH FILES AND PLOT w(t) and v(t) ----
# ####
# normalize_it <- function(x) x/(max(x))
# pdf(file = "./figures/all_wz_vz.pdf", width = 11, height = 11)
# par(mfrow=c(4,4),tcl=-0.2,mgp=c(2,0.5,0),mar=c(3,4,1,1),cex.lab=1.2, las=1)
# all_files <- list.files(paste0(results_dir))
# counter <- 1
# for(do_file in all_files){
#   file_info <- unlist(str_split(do_file,"_"))
#   state     <- file_info[1]
#   species   <- file_info[2]
#   
#   ##  Rename states with formatting for figures
#   if(state == "arizona")   {state <- "Arizona"}
#   if(state == "kansas")    {state <- "Kansas"}
#   if(state == "montana")   {state <- "Montana"}
#   if(state == "newmexico") {state <- "New Mexico"}
#   if(state == "idaho")     {state <- "Idaho"}
#   
#   tmp_wt <- readRDS(paste0(results_dir,do_file))[["wt"]]
#   tmp_vt <- readRDS(paste0(results_dir,do_file))[["vt"]]
#   tmp_u  <- readRDS(paste0(results_dir,do_file))[["size_bins"]]
#   
#   if(ncol(tmp_wt) != length(tmp_u)){
#     udiff <- tmp_u[2]-tmp_u[1]
#     tmp_u <- c((tmp_u[1] - udiff), tmp_u)
#   }
#   
#   mean_wt <- normalize_it(colMeans(tmp_wt[201:nrow(tmp_wt),]))
#   mean_vt <- normalize_it(colMeans(tmp_vt[201:nrow(tmp_vt),]))
#   
#   plot(tmp_u, mean_wt, type = "l", ylab = "Stable distributions", xlab = "Size, z")
#   lines(tmp_u, mean_vt, lty = 2)
#   title(paste(state, species, sep=", "))
#   if(counter == 1){
#     legend(-1, 0.99, legend = c(expression(italic("w(z)")),expression(italic("v(z)"))), lty = c(1,2), bty = "n", cex = 1.3)
#   }
# 
#   
#   counter <- counter+1
# }
# dev.off()
