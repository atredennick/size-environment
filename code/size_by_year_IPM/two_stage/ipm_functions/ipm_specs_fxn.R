##  ipm_param_fxn.R: Function for defining IPM specifications.

#' Defines IPM specifications.
#' @author Andrew Tredennick, Britta Teller, Peter Adler, Stephen Ellner
#' @param grow_data Data frame of growth data for calculating size ranges.
#' @return A list with IPM specifications, prenamed.
get_ipm_specs <- function(grow_data){
  ##  Simulation specs
  A       <- 10000 # area of 100cm x 100cm quadrat
  tlimit  <- 2200  # number of years to simulate
  burn_in <- 200   # years to cut before calculations
  
  ##  Size and grid specs
  rec_area <- 0.25 # WARNING! not log transformed! 
  bigM     <- 100 # c(75,75,50,50) # set matrix dimension for each species
  maxSize  <- round(max(grow_data$area.t0)) # in cm^2
  #maxSize <- 3000  # in cm^2: PSSP=225 HECO=202  POSE=260  ARTR=3000  
  minSize=0.2  #cm^2
  
  ##  Upper and lower sizes for integral
  # go WAY beyond the range of observed sizes, just to be safe 
  L <- log(0.25) - 1
  U <- log(maxSize) + 0.5   
  
  ##  Mesh specs
  # boundary points b and mesh points y. 
  b <- L+c(0:bigM)*(U-L)/bigM # chops up the size interval (L-U) into bigM-equal-sized portions.
  u <- 0.5*(b[1:bigM]+b[2:(bigM+1)]) # middle of each n-equal-sized portion
  h <- u[2]-u[1] # step size for midpoint rule. (see equations 4 and 5 in Ellner and Rees (2006) Am Nat.)
  
  ##  Variables for Wr (crowding/competition) approximation        
  b.r <- sqrt(exp(b)/pi)
  u.r <- sqrt(exp(u)/pi)
  
  tmp        <- range(u.r)
  size_range <- seq(tmp[1],tmp[2],length=50) # range across all possible sizes
  return(list(A          = A,
              tlimit     = tlimit,
              burn_in    = burn_in,
              rec_area   = rec_area,
              bigM       = bigM,
              maxSize    = maxSize,
              L          = L,
              U          = U,
              b          = b,
              u          = u,
              h          = h,
              b.r        = b.r,
              u.r        = u.r,
              size_range = size_range))
}

