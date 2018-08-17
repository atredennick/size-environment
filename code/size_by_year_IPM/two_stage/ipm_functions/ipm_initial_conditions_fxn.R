#' Function to generate initial population density vector
#' @author Britta Teller, Andrew Tredennick, Stephen Ellner
#' @param tlimit Number of iterations for IPM.
#' @param u Population vector.
#' @param surv_data Survival data frame to sample years.
#' @return A list of initial conditions, all named.
get_initial_conditions <- function(tlimit, u, h, A, surv_data){
  ##  Generate 10 seedlings 
  nt1 <- dnorm(u,mean=log(0.25),sd=0.001)
  nt1 <- 10 * nt1/sum(h*nt1)
  
  ##  Generate 3 adults
  nt2 <- dnorm(u,mean=5,sd=1) 
  nt2 <- 3 * nt2/sum(h*nt2)
  
  ##  Set up vector to record percent cover
  cov_save_kids      <- rep(NA,tlimit)
  cov_save_adults    <- rep(NA,tlimit)
  cov_save_kids[1]   <- sumCover(u,nt1,h,A)
  cov_save_adults[1] <- sumCover(u,nt2,h,A)
  
  ##  Set up matrix to store size distributions
  size_save_kids       <- matrix(NA,length(u),tlimit)
  size_save_adults     <- matrix(NA,length(u),tlimit)
  total_size           <- sum(nt1,nt2)
  size_save_kids[,1]   <- nt1/total_size
  size_save_adults[,1] <- nt2/total_size
  
  ##  Set up vector for initial densities 
  N_save_adults    <- rep(NA,tlimit)
  N_save_kids      <- rep(NA,tlimit)
  N_save_adults[1] <- sumN(nt1,h)
  N_save_kids[1]   <- sumN(nt2,h)
  
  ##  Sample years
  all_years <- as.numeric(as.character(unique(surv_data$year)))
  year_save <- sample(all_years, tlimit, replace = TRUE)
  return(
    list(
      nt1 = nt1,
      nt2 = nt2,
      cov_save_kids = cov_save_kids,
      cov_save_adults = cov_save_adults,
      size_save_kids = size_save_kids,
      size_save_adults = size_save_adults,
      total_size = total_size,
      N_save_kids = N_save_kids,
      N_save_adults = N_save_adults,
      year_save = year_save
    )
  )
}
