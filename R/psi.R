#' Anamorphosis function
#'
#' A function that computes the anamorphosis of four meta-Gaussian models (see Details)
#'
#' @param theta model parameters
#' @param x value to apply psi (Gaussian)
#' @param ym minimal value that can be observed
#' @param name name of the anamorphosis, one of \code{'power'},\code{'power-exp'}, \code{'quadratic-power'}, \code{'gp'}
#'
#' @return Precipitation amount from transforming x
#' @export
#'
#' @examples psi(c(-1,.2,.6,.1), 0.87, name='gp')
psi = function(theta, x, ym=0, name="power-exp"){
  switch(name,
         "power-exp" = ym+theta[3]*(exp(theta[2]*x^theta[4])-1),
         "power" = ym+theta[2]*x^theta[3],
         "quadratic-power"= ym+theta[2]*x^theta[4]+theta[3]*x^(2*theta[4]),
         "gp" = sapply(x, function(x) if (x<sqrt(abs(1/min(theta[3]*theta[4],0)))){ym+theta[2]*x^(1/theta[3])*exp(0.5*theta[4]*x^2)}else{NA})
  )
}
