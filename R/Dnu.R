#' D
#'
#' Tool function
#'
#' @param z input
#' @param nu input
#'
#' @importFrom fAsianOptions whittakerW
#' @export
Dnu=function(z,nu){2^(nu/2)*exp(-z^2/4)*(sqrt(pi)/gamma((1-nu)/2)*kummerM(z^2/2, -nu/2,1/2)-
                                           sqrt(2*pi)*z/gamma(-nu/2)*kummerM(z^2/2, (1-nu)/2,3/2))}
