################################################################################
##  vital_rate_ipm_fxns.R: This script contains functions for estimating
##  growth, survival, and recruitment based on fitted vital rate functions and
##  the current state of the population in the Integral Projection Model (IPM).
##
##
##  ____________________________________________________________________________
##  Author: Andrew Tredennick (atredenn@gmail.com)
##  Contributors: Peter Adler, Britta Teller, Stephen Ellner
##  Date created: August 28, 2017
################################################################################


#' Predict growth based on crowding
#' 
#' @author Andrew Tredennick, Britta Teller, Peter Adler, Stephen Ellner
#' @param u Size class vector.
#' @param v Size class vector, again.
#' @param W Crowding vector for intraspecific density dependence.
#' @param Gpars Growth regression parameters.
#' @param do_year The year we are working with.
#' @return Growth kernel.
G <- function(u,v,W,Gpars,do_year,grow_allsize_years){
  Gpars <- subset(Gpars, year_id == do_year)
  sxy <- grow_allsize_years[which(as.numeric(row.names(grow_allsize_years)) == do_year),]
  mu    <- Gpars$intercept + (Gpars$W)*W + Gpars$logarea.t0*u + (Gpars$logarea.W)*(u*W) + sxy
  # sigma <- exp(Gpars$sigma.a+Gpars$sigma.b*mu)
  # return( dT(x=v, mu=mu, sigma=sigma, nu=exp(Gpars$sigma.nu)) )
  
  sigma2 <- Gpars$sigma.a*exp(Gpars$sigma.b*mu)

  out <- dnorm(v,mu,sqrt(sigma2))
  return(out)
}



#' Predict survival based on crowding
#' 
#' @author Andrew Tredennick, Britta Teller, Peter Adler, Stephen Ellner
#' @param u Size class vector.
#' @param v Size class vector, again.
#' @param W Crowding vector for intraspecific density dependence.
#' @param Spars Survival regression parameters.
#' @param do_year The year we are working with.
#' @return Survival kernel.
S <- function(u,v,W,Spars,do_year,surv_allsize_years){
  Spars <- subset(Spars, year_id == do_year)
  sxy <- surv_allsize_years[which(as.numeric(row.names(surv_allsize_years)) == do_year),]
  mu <- Spars$intercept + (Spars$W)*W + Spars$logarea*u + (Spars$logarea.W)*(u*W) + sxy

  return(inv.logit(mu))
}



#' Import empirical size distribution of new recruits
#' @author Andrew Tredennick
#' @param data_dir Directory with data.
#' @param state The state we are dealing with.
#' @param species The species we are dealing with.
#' @return A vector (length = 2) with recSize and recVariance.
get_rec_size <- function(data_dir, state, species){
  infile <- paste0(data_dir,state,"/speciesData/",species,"/recSize.csv")
  if(file.exists(infile)==FALSE){
    infile <- paste0(data_dir,state,"/speciesData/",species,"/recSize_reduced.csv")
  }
  if(file.exists(infile)==FALSE){
    stop("Recruit size data file not found. Check filename.")
  }
  rec_size <- read.csv(infile)
  return(c(sizeMean = mean(log(rec_size$area)), 
           sizeVar  = var(log(rec_size$area))))
}



#' Predict recruitment based on total cover
#' 
#' @author Andrew Tredennick, Britta Teller, Peter Adler, Stephen Ellner
#' @param params Recruitment regression parameters.
#' @param cover Cover vector, in m^2 per m^2.
#' @param A Area being modeled.
#' @param do_year The year we are working with.
#' @return Vector of number of recruits per area produced.
get_recs_per_area <- function(params,cover,A,do_year){
  Rpars <- subset(params, year == do_year)
  ##  Cover is in m^2 per m^2; convert to % scale:
  cover2 <- cover*100
  
  ##  Calculate recruits
  mu <- cover2*exp(Rpars$intcpt.mu + Rpars$intcpt.yr + sqrt(cover2)%*%Rpars$dd) 
  if(sum(is.na(mu))>0) browser() # stop for errors
  rpa <- mu/(cover*A)  # convert from number recruits to recruits per cm^2
  return(rpa)
}

#' Fecundity function, expected number of recruits of size y produced by a size x individual
#' The size distribution of recruits is on the log scale
#' @author Andrew Tredennick, Britta Teller, Peter Adler, Stephen Ellner
#' @param u Current population vector.
#' @param params Empirical recruitment parameters: mean size and size variance.
#' @param rpa Recruits per area, from get_recs_per_area.
f <- function(u,v,Rpars,rpa){ 
  nRecruits <-  rpa*exp(u)
  
  # recSizeDist <- dnorm(u,mean=log(0.25),sd=0.2); 
  # recSizeDist <- 0.25*recSizeDist/sum(h*exp(u)*recSizeDist); 
  #sum(h*exp(u)*recSizeDist)/recArea; # needs to equal 1 
  
  ##peter's older distribution that isn't long-tailed enough
  tmp=dnorm(v,Rpars[1],sqrt(Rpars[2]))/(1-pnorm(-1.61,Rpars[1],sqrt(Rpars[2])))
  #number recruits of each size 
  
  f=nRecruits*tmp
  
  ##  Probability of producing a seedling of size v
  #tmp <- dnorm(v,Rpars["sizeMean"],sqrt(Rpars["sizeVar"]))/(1-pnorm(-1.61,Rpars["sizeMean"],sqrt(Rpars["sizeVar"])))
  ##  Number recruits of each size 
  return(f)
}   
