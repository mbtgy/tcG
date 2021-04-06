#' QQplot for the extended GP model
#'
#' Called by \code{res.plot}, we don't advise you to use it.

#' @param y vector of data
#' @param qfit list obtained with \code{fit.extGP} (fitted quantiles)
#' @param q.L list obtained with \code{fit.extGP} (lower intervall of the bootstrap quantiles)
#' @param q.U list obtained with \code{fit.extGP} (upper intervall of the bootstrap quantiles)
#' @param cols vector of colors for the barplot
#' @param zoomx vector of length 2, range to zoom on the x axis
#' @param zoomy vector of length 2, range to zoom on the y axis
#' @param legend boolean, should there be a legend?
#' @param main title of the plot
#' @param add boolean, should the plot be added to the current plot?
#'
#' @importFrom grDevices rgb
#' @importFrom graphics polygon lines legend abline
#'
#' @export
#'
extGP.qqplot = function(y, qfit, q.L, q.U, cols=c(rgb(.6,.3,0),rgb(.6,.3,0,.5), rgb(.6,.3,0,.1)),
                        zoomx=range(y,na.rm=T), zoomy=range(y,na.rm=T), legend=TRUE,
                        main="QQ plot of positive rainfall", add=FALSE){
  if (!add) plot(0, 0, asp=1, xlab="Empirical quantiles",ylab="Fitted quantiles", main=main,
                 ylim=zoomy, xlim=zoomx, type="n")
  if (length(q.L)!=0){
    M=cbind(c(sort(y), sort(y,T)), c(q.L, q.U[length(q.U):1]))
    M=M[!duplicated(M),]
    polygon(x=M[,1], col=cols[3], border=cols[2],
            y=M[,2], lty=2)
  }
  M=cbind(sort(y), qfit)
  M=M[!duplicated(M),]
  lines(M[,1], M[,2], lty=2, type="b", pch=20, col=cols[1])
  abline(0,1)
  if (legend) legend("bottomright", "extGP", text.col=cols[1], bty="n")
}
