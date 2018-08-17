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



# Adult growth ------------------------------------------------------------

#' Predict growth based on current population state
#' 
#' @author Andrew Tredennick, Britta Teller, Peter Adler, Stephen Ellner
#' @param u Population size distribution vector.
#' @param v Population size distribution vector (copy of u).
#' @param W Crowding vector for intraspecific density dependence.
#' @param Gpars Growth regression parameters.
#' @param do_year The year we are working with.
#' @return Growth kernel.
G <- function(u,v,W,Gpars,do_year){
  Gpars <- subset(Gpars, year_id == do_year)
  
  mu <- Gpars$intercept +     # intercept
    (Gpars$W)*W +           # crowding effect
    Gpars$logarea.t0*u +      # size effect
    (Gpars$logarea.W)*(u*W) # size*crowding effect
  
  ##  PBA Normal distribution version
  sigma2 <- Gpars$sigma.a*exp(Gpars$sigma.b*mu)
  return(dnorm(v,mu,sqrt(sigma2)))
  
  ## SPE T-Distribution version
  # sigma <- exp(Gpars$sigma.a+Gpars$sigma.b*mu)
  # return(dT(x=v,mu=mu,sigma=sigma, nu=exp(Gpars$sigma.nu)))
}



# Seedling growth ---------------------------------------------------------

#' Predict growth of seedlings
#' 
#' @author Andrew Tredennick, Britta Teller, Peter Adler, Stephen Ellner
#' @param u Population size distribution vector.
#' @param v Population size distribution vector (copy of u).
#' @param pars Probability of growing and empirical mean and sd of seedling growth.
#' @param kidsize The log size of seedlings, defined by user/data.
#' @return Size distribution vector, length of u 
get_growth_kids <- function(u,v,pars,kidsize) {
  kidbin        <- 1 + sum(v<kidsize)
  h             <- v[2]-v[1] 
  part2         <- dnorm(v,mean=pars[2],sd=pars[3])
  part2         <- (1-pars[1])*part2/sum(h*part2)
  part1         <- rep(0,length(v))
  part1[kidbin] <- pars[1]/h
  return(part1+part2)
}
# get_growth_kids <- function(u,v,pars){
#   return( pars[1]*dnorm(v,mean=log(0.25),sd=0.01)+(1-pars[1])*dnorm(v,mean=pars[2],sd=pars[3]) )
# }




# Adult survival ----------------------------------------------------------

#' Predict survival based on current population state
#' 
#' @author Andrew Tredennick, Britta Teller, Peter Adler, Stephen Ellner
#' @param u Population size distribution vector.
#' @param v Population size distribution vector (copy of u).
#' @param W Crowding vector for intraspecific density dependence.
#' @param Spars Survival regression parameters.
#' @param do_year The year we are working with.
#' @return Growth kernel.
S <- function(u,v,W,Spars,do_year){
  Spars <- subset(Spars, year_id == do_year)
  
  mu <- Spars$intercept +     # intercept
    (Spars$W)*W +           # crowding effect
    Spars$logarea*u +         # size effect
    (Spars$logarea.W)*(u*W) # size*crowding effect
  
  return(inv.logit(mu))
}




# Seedling survival -------------------------------------------------------

#' Predict survival of seedlings.
#' 
#' @author Andrew Tredennick, Britta Teller, Peter Adler, Stephen Ellner
#' @param u Population size distribution vector.
#' @param v Population size distribution vector (copy of u).
#' @param W Crowding vector for intraspecific density dependence.
#' @param Spars Survival regression parameters for seedlings.
#' @param do_year The year we are working with.
#' @return Growth kernel.
get_survival_kids <- function(u,v,W,Spars,do_year){
  Spars <- subset(Spars, year_id == do_year)
  
  mu <- Spars$intercept + # intercept
    (Spars$W)*W         # crowding effect
  
  return(inv.logit(mu))
}




# New recruit size --------------------------------------------------------

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




# Number of new recruits --------------------------------------------------

#' Predict recruitment based on current population state
#' 
#' @author Andrew Tredennick, Britta Teller, Peter Adler, Stephen Ellner
#' @param params Recruitment regression parameters.
#' @param cover Cover vector, in m^2 per m^2.
#' @param do_year The year we are working with.
#' @return Vector of number of recruits.
get_rec_number <- function(params,cover,do_year){
  Rpars <- subset(params, year == do_year)
  ##  Cover is in m^2 per m^2; convert to % scale:
  cover2 <- cover*100
  
  ##  Calculate recruits
  mu <- cover2*exp(Rpars$intcpt.mu + Rpars$intcpt.yr + sqrt(cover2)%*%Rpars$dd) 
  if(sum(is.na(mu))>0) browser() # stop for errors
  return(mu) # return number of recruits
}




# Recruitment density -----------------------------------------------------

#' Fecundity function, expected number of recruits of size y produced by a size x individual
#' The size distribution of recruits is on the log scale
#' @author Andrew Tredennick, Britta Teller, Peter Adler, Stephen Ellner
#' @param u Population size distribution vector.
#' @param v Population size distribution vector (copy of u).
#' @param Rpars Empirical recruitment parameters: mean size and size variance.
#' @param rpa Number of new recruits, from get_recs_per_area.
#' @return A scalar representing the density of new seedlings for the nt1 vector kidsize bin.
f <- function(u,v,Rpars,rpa){ 
  recSizeDist <- dnorm(log(0.25),mean=log(0.25),sd=0.01);
  recSizeDist <- recSizeDist/sum(h*recSizeDist);
  return(rpa*recSizeDist)
}   


#' Fecundity function, expected number of recruits of size y produced by a size x individual
#' The size distribution of recruits is on the log scale
#' @author Andrew Tredennick, Britta Teller, Peter Adler, Stephen Ellner
#' @param u Current population vector.
#' @param params Empirical recruitment parameters: mean size and size variance.
#' @param rpa Recruits per area, from get_recs_per_area.
f2 <- function(u,v,Rpars,rpa){ 
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

