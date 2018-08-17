##  Kernel and projection matrix functions, based on vital rate regressions.
##  Also includes some other functions for summing cover and abundances over
##  the population state vector.

#' Calculate survival*growth kernel values, called from make.K.matrix
#' @author Peter Adler, Stephen Ellner
#' @param v Vector of size classes
#' @param u Vector of size classes
#' @param muW Expanded vector of crowding (see expandW function)
#' @param Gpars Growth regression parameters, full dataframe
#' @param Spars Survival regression parameters, full dataframe
#' @param do_year The year we want to grab for random year effects
#' @return A u*v sized projection matrix
make.P.values <- function(v,u,muW,Gpars,Spars,do_year){
  S(u, v , muW, Spars, do_year)*G(u, v, muW, Gpars, do_year) 
}


#' Calculate fecundity kernel values, called from make.K.matrix
#' @author Peter Adler, Stephen Ellner
#' @param v Vector of size classes
#' @param u Vector of size classes
#' @param Rpars Empirical recruitment parameters for recruit size and variance
#' @param rpa Calculated recruits per area
#' @param do_year The year we want to grab for random year effects
#' @return A u*v sized projection matrix
make.F.values <- function(v,u,Rpars,rpa,do_year){
  f(u, v, Rpars, rpa)
}


#' Calculate combined kernel
#' @author Peter Adler, Stephen Ellner
#' @param v Vector of size classes
#' @param muW Vector of crowding
#' @param Rpars Empirical recruitment parameters for recruit size and variance
#' @param rpa Calculated recruits per area
#' @param Gpars Growth regression parameters, full dataframe
#' @param Spars Survival regression parameters, full dataframe
#' @param do_year The year we want to grab for random year effects
#' @return A u*v sized projection matrix
make.K.matrix <- function(v,muW,Rpars,rpa,Gpars,Spars,do_year,fecundity_cutoff) {
  muW <- expandW(v,muW)
  
  ##  Get fecundity kernel first
  F.matrix <- outer(v,v,make.F.values,Rpars,rpa,do_year)

  ##  Truncate fecundity to that of max observed size
  for(i in fecundity_cutoff:100){
    F.matrix[,i] <- F.matrix[,fecundity_cutoff]
  }
  
  ##  Now get survival*growth kernel
  P.matrix <- outer(v,v,make.P.values,muW,Gpars,Spars,do_year)
  
  ##  Combine for population growth kernel
  K.matrix <- F.matrix + P.matrix
  return(h*K.matrix)
}


#' Calculate density of a t-distribution
#' @author Stephen Ellner
#' @param x Vector over which to calculate density
#' @param mu Mean
#' @param sigma scale parameter ???
#' @param nu rate parameter ???
#' @return A vector of densities for each size class.
dT <- function(x,mu,sigma,nu) {
  a=1/2; b=nu/2; B = gamma(a)*gamma(b)/gamma(a+b); 
  fac1 = 1/(sigma*B*sqrt(nu)); 
  fac2 = 1 + ((x-mu)^2)/(nu*sigma^2);
  pwr = -(1+nu)/2; 
  return(fac1*(fac2^pwr))    
}



#' Function to sum total % cover of each species
#' @author Britta Teller, Peter Adler, Stephen Ellner
#' @param u Vector of size classes
#' @param nt Population abundance, binned.
#' @param h Midpoint.
#' @param A Area being modeled.
#' @return Summed cover over the area, cm^2 per cm^2
sumCover <- function(u,nt,h,A){
  return( h*sum(nt*exp(u))/A )
} 

#' Function to sum total N 
#' @author Britta Teller, Peter Adler, Stephen Ellner
#' @param nt Population abundance, binned.
#' @param h Midpoint.
#' @return A numeric scalar of total abundance.
sumN <- function(nt,h){
  return( h*sum(nt) )
}



#' Make a nice plot of the projection kernel
#' @author Stephen Ellner
#' @return A plot
matrix.image=function(A, x=NULL, y=NULL, col=rainbow(100,start=0.67,end=0),
             bw=FALSE, do.contour=FALSE, do.legend=TRUE,...) {
 if(do.legend) layout(mat=cbind(matrix(1,5,5),rep(2,5)));
 par(mar=c(6,5,3,2));
 if(is.null(x)) x=1:ncol(A);
 if(is.null(y)) y=1:nrow(A);
 nx=length(x); ny=length(y);
 x1=c(1.5*x[1]-0.5*x[2],1.5*x[nx]-0.5*x[nx-1]);
 y1=c(1.5*y[1]-0.5*y[2],1.5*y[ny]-0.5*y[ny-1]);
 if(bw) col=grey( (200:50)/200 );
 image(list(x=x,y=y,z=t(A)),xlim=x1,ylim=rev(y1),col=col,cex.axis=1.5,cex.lab=1.5,bty="u",...);
 abline(v=range(x1)); abline(h=range(y1));
 if(do.contour) contour(x,y,t(A),nlevels=5,labcex=1.2,add=TRUE);

 if(do.legend) {
    l.y=seq(min(A),max(A),length=100);
    par(mar=c(6,2,3,1))
    image(list(x=1:2,y=l.y,z=rbind(l.y,l.y)),col=col,bty="o",xaxt="n",yaxt="n");
    axis(side=2,cex.axis=1.5,at=pretty(seq(min(A),max(A),length=10)));
 }
}
