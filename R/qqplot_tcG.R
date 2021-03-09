#' QQplot for meta-Gaussian models
#'
#' Called by \code{res.plot}, we don't advise you to use it.
#'
#' @param y vector of data
#' @param names name(s) of the anamorphosis that were fitted, a subset of \code{c('gp', 'power', 'power-exp', 'quadratic-power')}
#' @param qfit list obtained with \code{fit.tcG} (fitted quantiles)
#' @param q.L list obtained with \code{fit.tcG} (lower intervall of the bootstrap quantiles)
#' @param q.U list obtained with \code{fit.tcG} (upper intervall of the bootstrap quantiles)
#' @param cols matrix of colors for the barplot
#' @param zoomx vector of length 2, range to zoom on the x axis
#' @param zoomy vector of length 2, range to zoom on the y axis
#' @param legend boolean, should there be a legend?
#' @param main title of the plot
#'
#' @importFrom grDevices rgb
#' @importFrom graphics polygon lines legend abline
#'
#' @export
#'

tcG.qqplot = function(y, names, qfit, q.L, q.U, cols=cbind(c(rgb(0,0,1), rgb(1,0,0), rgb(0,1,0), rgb(1,0,1)),
                                                           c(rgb(0,0,1,0.5), rgb(1,0,0,0.5), rgb(0,1,0,0.5), rgb(1,0,1,0.5)),
                                                           c(rgb(0,0,1,0.1), rgb(1,0,0,0.1), rgb(0,1,0,0.1), rgb(1,0,1,0.1))),
                      zoomy=range(y,na.rm=T),legend=TRUE, zoomx=range(y,na.rm=T), main="QQ plot"){
  plot(0, 0, asp=1, xlab="Empirical quantiles",ylab="Fitted quantiles", main=main,
       ylim=zoomy, xlim=zoomx, type="n",cex.main=1)
  i=1
  for (nm in names){
    if (length(q.L)!=0){
      M=cbind(c(sort(y), sort(y,T)), c(q.L[[nm]], q.U[[nm]][length(q.U[[nm]]):1]))
      M=M[!duplicated(M),]
      polygon(x=M[,1], col=cols[i,3], border=cols[i,2],
              y=M[,2], lty=2)
    }
    M=cbind(sort(y), qfit[[nm]])
    M=M[!duplicated(M),]
    lines(M[,1], M[,2], lty=2, type="b", pch=20, col=cols[i,1])
    i=i+1
  }
  abline(0,1)
  if (legend) legend("topleft", names, text.col=cols[1:(i-1),1], bty="n")
}
