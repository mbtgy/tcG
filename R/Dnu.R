#' D
#'
#' Tool function
#'
#' @param z input
#' @param nu input
#'
#' @importFrom fAsianOptions whittakerW
#' @export
Dnu=function(z,nu){2^(1/4+nu/2)*z^(-1/2)*whittakerW(z^2/2, kappa=1/4+nu/2,mu=-1/4)}
