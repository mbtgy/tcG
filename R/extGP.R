#' Extended Generalized Pareto distribution
#'
#' Density, distribution function, quantile function and random generation of an extended Generalised Pareto distribution with parameters theta.
#'
#' @param y vector of quantiles (rainfall)
#' @param theta model parameters (sigma, alpha, xi)
#' @param ym minimal value that can be observed
#'
#' @importFrom evd dgpd pgpd
#' @return \code{dextGP} gives the density, \code{pextGP} gives the distribution function, \code{qextGP} gives the quantile function and \code{rextGP} generates random deviates.
#' @export

pextGP=function(y,theta, ym=0){
  p=(1-pmax((1+theta[3]*(y-ym)/theta[1]),0)^(-1/theta[3]))^(theta[2])  #renvoie NaN si x<theta[4]
  p[is.na(p)]=0
  return(p)
}

#' @param p vector of probabilities
#' @param step discretization step
#' @describeIn pextGP quantile function
#' @export
qextGP=function(p,theta,ym=0, step=0){
  qcont=ym+theta[1]/theta[3]*((1-p^(1/theta[2]))^(-theta[3])-1)  #sans discrétisation
  if (step!=0) qcont=step*floor(qcont/step)
  return(qcont)
}

#' @describeIn pextGP density function
#' @export
dextGP=function(y,theta,ym=0,step=0){
  if (step==0){
  kappa=theta[2];sigma=theta[1];xi=theta[3]
  return(kappa * pgpd(y+ym, scale = sigma, shape = xi)^(kappa - 1)*
           dgpd(y+ym, scale = sigma, shape = xi))
  }else{
    pextGP(y+step, theta, ym)-pextGP(y, theta, ym)
  }
}
