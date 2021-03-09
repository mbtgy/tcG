#' Moments of GP meta-Gaussian model
#'
#' Computes the moments of the GP meta-Gaussian model, all distribution (momtcG) and positive part (momtcG.pos)
#'
#' @param p order of the moment
#' @param mu mean of the Gaussian
#' @param sigma scale parameter
#' @param alpha lower tail parameter
#' @param xi upper tail parameter
#' @param ym minimal value that can be observed
#'
#' @importFrom stats pnorm
#'
#' @return Moment of order p of the GP meta-Gaussian distribution
#' @export
#'
momtcG = function(p, mu, sigma, alpha, xi, ym=0){
  (1-pnorm(-mu))*momtcG.pos(p, mu, sigma, alpha, xi)+(-ym)^p*pnorm(-mu)
}



#' @export
#' @describeIn momtcG Moment of the positive part of the GP meta-Gaussian distribution
momtcG.pos = function(p, mu, sigma, alpha, xi){
  sigma^p*(1-xi*p)^(-(p/alpha+1)/2)/(sqrt(2*pi)*(1-pnorm(-mu)))*
    exp(mu^2/2*(1/(2*(1-xi*p))-1))*
    gamma(p/alpha+1)*Dnu(z=-mu/sqrt(1-xi*p), nu=-(p/alpha+1))
}

