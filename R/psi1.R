#' Inverse of the anamorphosis
#'
#' A function that computes the inverse of the anamorphosis of four meta-Gaussian models (see Details)
#'
#' @param theta model paramters
#' @param y value to apply psi^{-1} (rainfall)
#' @param ym minimal value that can be observed
#' @param name name of the anamorphosis, one of \code{'power'},\code{'power-exp'}, \code{'quadratic-power'}, \code{'gp'}
#'
#' @return Returns the inverse of psi
#' @importFrom VGAM lambertW
#' @export
#'
#' @examples psi1(c(-1,.2,.6), 0.1, name='power')
psi1 = function(theta, y, ym=0, name="power-exp"){
  switch(name,
         "power-exp" = (log((y-ym)/theta[3]+1)/theta[2])^(1/theta[4]),
         "power" = ((y-ym)/theta[2])^(1/theta[3]),
         "quadratic-power"= ((sqrt(theta[2]^2+4*theta[3]*(y-ym))-theta[2])/(2*theta[3]))^(1/theta[4]),
         "gp" = sqrt(lambertW(theta[4]*theta[3]*((y-ym)/theta[2])^(2*theta[3]))/(theta[3]*theta[4]))
  )
}



