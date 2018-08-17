##  ipm_crowding_fxns.R: Functions for estimating crowding effects for the IPM.

#' Get variables from data to pass to approximation function.
#' @author Andrew Tredennick, Britta Teller, Peter Adler, Stephen Ellner
#' @param data_dir Directory where dataframe of distance weights calculated from data can be found.
#' @param state The state we are dealing with.
#' @param species The species we are dealing with.
#' @param b.r
#' @return A list of variables necessary to approximate crowding in IPM.
get_crowd_vars <- function(data_dir, state, species, b.r){
  ##  Get empirical distance weights for state
  dist_wts <- read.csv(paste0(data_dir,state,"/speciesdata/",state,"DistanceWeights.csv"))
  
  ##  Make dist_wts constant over first annulus (0 to 2 cm)
  dist_wts <- rbind(dist_wts,c(1,1,1,1,0))
  dist_wts <- dist_wts[order(dist_wts$midRings),]
  Wfuns <- approxfun(x=dist_wts$midRings, y=dist_wts[,species], rule=2) 
  
  tmp <- which(dist_wts[,species]==0)
  if(length(tmp)>0){
    zc <- dist_wts$midRings[min(tmp)]
  }else{
    zc <- 150  # max distance
  }
  
  return(list(b.r=b.r,Wfuns=Wfuns,zc=zc,h=h))
}



#' Make a crowding matrix for the IPM.
#' @author Andrew Tredennick, Britta Teller, Peter Adler, Stephen Ellner
#' @param Wvars List of crowding variables from get_crowd_vars.
#' @param u Vector of size clases.
#' @param A Size of simulation area.
makeW <- function(Wvars, u, A, Cr){
  Wfuns <- Wvars$Wfuns
  
  wrii <- function(r) {
    lower = r; upper=Wvars$zc; 
    pts= seq(lower,upper,length=501); dx=pts[2]-pts[1]; 
    vals = pts*Wfuns(pts)*Cr(pts-r); 
    int = dx*sum(vals[2:500])+0.5*dx*(vals[1]+vals[501]) 
    return(2*pi*int)
  }   
  
  u.r  <- sqrt(exp(u)/pi) # radius of genet of each log size u
  nu = length(u.r); Wmat = numeric(nu); 
  for(j in 1:nu){Wmat[j] <- wrii(u.r[j])/A}
  
  return(Wmat)
}

# makeW <- function(Wvars, u, A, Cr){
#   Wfuns <- Wvars$Wfuns
# 
#   # if W is continuous
#   wrii <- function(r){
#     return(2*pi*integrate(function(z) z*Wfuns(z)*Cr(z-r), r, Wvars$zc)$value)
#   }
#   Wrii <- Vectorize(wrii,vectorize.args="r")
# 
#   u.r  <- sqrt(exp(u)/pi) # radius of genet of each log size u
#   Wmat <- Wrii(u.r)/A
# 
#   return(Wmat)
# }


#' Expand the crowding vector for outer product usage
#' @author Stephen Ellner
#' @param u Vector of size classes
#' @param W Vector of crowding
#' @return A long vector of size W*u for all combinations of W and u
expandW <- function(u, W){
  if(length(W)!=length(u)) stop("Check length of W")
  W <- rep(W, each=length(u))
  return(W)
}



#' Calculate size-dependent crowding, assuming no overlap
#' @author Stephen Ellner
#' @param r Radius of plant of a given area
# wrii <- function(r) {
#   return(2*pi*integrate(function(z) z*Wvars$z*Cr(z-r),r,Wvars$zc)$value)
# }
# Wrii <- Vectorize(wrii,vectorize.args="r")

