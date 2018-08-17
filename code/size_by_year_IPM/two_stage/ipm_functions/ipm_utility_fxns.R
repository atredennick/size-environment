##  Kernel and projection matrix functions, based on vital rate regressions.
##  Also includes some other functions for summing cover and abundances over
##  the population state vector.

dT <- function(x,mu,sigma,nu) {
  a=1/2; b=nu/2; B = gamma(a)*gamma(b)/gamma(a+b); 
  fac1 = 1/(sigma*B*sqrt(nu)); 
  fac2 = 1 + ((x-mu)^2)/(nu*sigma^2);
  pwr = -(1+nu)/2; 
  return(fac1*(fac2^pwr))    
}



## combined kernel
make.P.values=function(v,u,muW, #state variables
                       Gpars,Spars,do_year,G,S)  #growth arguments
{
  S(u, v , muW, Spars, do_year)*G(u, v, muW, Gpars, do_year) 
}


# Function to make iteration matrix based only on mean crowding
make.P.matrix=function(v,Wvec,Gpars,Spars,do_year=NULL,G,S) {
  muW <- expandW(v,Wvec)
  P.matrix <- outer(v,v,make.P.values,muW,Gpars,Spars,do_year,G,S)
  return(h*P.matrix)
}





#' Function to sum total % cover of each species
#' @author Britta Teller, Peter Adler, Stephen Ellner
#' @param u Population state vector.
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
sumN <- function(nt,h){
  return( h*sum(nt) )
}


####
#### UNUSED FUNCTIONS
####

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
