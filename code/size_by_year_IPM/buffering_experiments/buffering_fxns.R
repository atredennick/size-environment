################################################################################
##  buffering_fxns.R: This script contains functions to induce size*year
##  buffering. Takes fitted parameter values from a model with only a year
##  effect on the intercept and creates scenarios of size*year interactions.
##  The three scenarios are:
##    1. Small plants more variable than large plants
##    2. Small plants less variable than large plants
##    3. Small and large plants have opposing year effects
##
##  Examples of these can be seen by running the BufferingExperiments.R script
##
##  Authors: Andrew Tredennick (atredenn@gmail.com)
##           Stephen Ellner
##  Date created: 11-6-2017
################################################################################

make_small_variable <- function(sizes, mod_fit, view_plot = FALSE){
  ##  Baseline MS anomaly 
  alpha  <- ranef(mod_fit)$fyear[,1]
  Ab     <-  sqrt(mean(alpha^2))
  Nyears <- length(alpha)
  
  ##  Parameters for defining SxC models 
  z    <- sizes
  zL   <- quantile(z,0.05)
  zU   <- quantile(z,0.95)
  c    <- 3/(zU-zL); 
  zbar <- ztilde <- (zL + zU)/2
  
  ##  Compute anomalies for q=1
  anomalies <- numeric(0) 
  for(j in 1:Nyears) {
    aj <-  alpha[j]*exp(-c*(z-ztilde)) 
    anomalies <- c(anomalies,aj)
  }   
  A <- sqrt(mean(anomalies^2))
  q <- Ab/A # scale up to match baseline RMS anomaly. 
  
  ##  Create Size*Year effect matrix
  size_year <- matrix(NA,Nyears,length(z))
  for(j in 1:Nyears){
    size_year[j,] <- q*alpha[j]*exp(-c*(z-ztilde))
  }
  
  if(view_plot == TRUE){
    ##  Plot check
    a0 <- fixef(mod_fit)[1]
    a1 <- fixef(mod_fit)[2]
    a2 <- fixef(mod_fit)[3]
    a3 <- fixef(mod_fit)[4]
    wbar <- mean(gdata$W)
    yhat <- matrix(NA,length(z),Nyears)
    for(j in 1:Nyears){
      yhat[,j] <- a0 + a1*z + a2*wbar + a3*wbar*z + q*alpha[j]*exp(-c*(z-ztilde)) 
    }
    matplot(z,yhat,type="l",lty=2,col="black",main="Small are more variable",xlab="Size z"); 
    points(z,a0 + a1*z + a2*wbar + a3*wbar*z,type="l",lty=1,col="blue",lwd=2); 
  }
  row.names(size_year) <- as.numeric(row.names(coef(mod_fit)$fyear))
  return(size_year)
}


make_large_variable <- function(sizes, mod_fit, view_plot = FALSE){
  ##  Baseline MS anomaly 
  alpha  <- ranef(mod_fit)$fyear[,1]
  Ab     <-  sqrt(mean(alpha^2))
  Nyears <- length(alpha)
  
  ##  Parameters for defining SxC models 
  z    <- sizes
  zL   <- quantile(z,0.05)
  zU   <- quantile(z,0.95)
  c    <- 3/(zU-zL); 
  zbar <- ztilde <- (zL + zU)/2
  
  ##  Compute anomalies for q=1
  anomalies <- numeric(0) 
  for(j in 1:Nyears) {
    aj <-  alpha[j]*exp(c*(z-ztilde)) 
    anomalies <- c(anomalies,aj)
  }   
  A <- sqrt(mean(anomalies^2))
  q <- Ab/A # scale up to match baseline RMS anomaly. 
  
  ##  Create Size*Year effect matrix
  size_year <- matrix(NA,Nyears,length(z))
  for(j in 1:Nyears){
    size_year[j,] <- q*alpha[j]*exp(c*(z-ztilde))
  }
  
  if(view_plot == TRUE){
    ##  Plot check
    a0 <- fixef(mod_fit)[1]
    a1 <- fixef(mod_fit)[2]
    a2 <- fixef(mod_fit)[3]
    a3 <- fixef(mod_fit)[4]
    wbar <- mean(gdata$W)
    yhat <- matrix(NA,length(z),Nyears)
    for(j in 1:Nyears){
      yhat[,j] <- a0 + a1*z + a2*wbar + a3*wbar*z + q*alpha[j]*exp(c*(z-ztilde)) 
    }
    matplot(z,yhat,type="l",lty=2,col="black",main="Large are more variable",xlab="Size z"); 
    points(z,a0 + a1*z + a2*wbar + a3*wbar*z,type="l",lty=1,col="blue",lwd=2); 
  }
  
  row.names(size_year) <- as.numeric(row.names(coef(mod_fit)$fyear))
  return(size_year)
}


make_sign_flip <- function(sizes, mod_fit, view_plot = FALSE){
  ##  Baseline MS anomaly 
  alpha  <- ranef(mod_fit)$fyear[,1]
  Ab     <-  sqrt(mean(alpha^2))
  Nyears <- length(alpha)
  
  ##  Parameters for defining SxC models 
  z    <- sizes
  zL   <- quantile(z,0.05)
  zU   <- quantile(z,0.95)
  c    <- 3/(zU-zL); 
  zbar <- ztilde <- (zL + zU)/2
  dz   <- (zU-zL)
  
  ##  Compute anomalies for q=1
  anomalies <- numeric(0) 
  for(j in 1:Nyears) {
    aj <-  2*alpha[j]*(z-zbar)/dz 
    anomalies <- c(anomalies,aj)
  }   
  A <- sqrt(mean(anomalies^2))
  q <- Ab/A # scale up to match baseline RMS anomaly. 
  
  ##  Create Size*Year effect matrix
  size_year <- matrix(NA,Nyears,length(z))
  for(j in 1:Nyears){
    size_year[j,] <- q*2*alpha[j]*(z-zbar)/dz 
  }
  
  if(view_plot == TRUE){
    ##  Plot check
    a0 <- fixef(mod_fit)[1]
    a1 <- fixef(mod_fit)[2]
    a2 <- fixef(mod_fit)[3]
    a3 <- fixef(mod_fit)[4]
    wbar <- mean(gdata$W)
    yhat <- matrix(NA,length(z),Nyears)
    for(j in 1:Nyears){
      yhat[,j] <- a0 + a1*z + a2*wbar + a3*wbar*z + q*2*alpha[j]*(z-zbar)/dz 
    }
    matplot(z,yhat,type="l",lty=2,col="black",main="Large are more variable",xlab="Size z"); 
    points(z,a0 + a1*z + a2*wbar + a3*wbar*z,type="l",lty=1,col="blue",lwd=2); 
  }
  
  row.names(size_year) <- as.numeric(row.names(coef(mod_fit)$fyear))
  return(size_year)
}