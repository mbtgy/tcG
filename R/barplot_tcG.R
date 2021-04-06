#' Barplot for meta-Gaussian fitted model
#'
#' Called by \code{res.plot}, we don't advise you to use it.
#'
#' @param y vector of data
#' @param theta vector of model parameters
#' @param name name of the anamorphosis, one of \code{'power'},\code{'power-exp'}, \code{'quadratic-power'}, \code{'gp'}
#' @param ym minimal value that can be observed
#' @param step discretization step
#' @param cols vector of length 2, colors for the barplot
#' @param zoom vector of length 2, range to zoom on the x axis
#' @param legend boolean, should there be a legend?
#' @param main title of the plot
#'
#' @importFrom grDevices rgb
#' @importFrom graphics barplot legend
#'
#' @export
#'
tcG.barplot = function(y, theta, name="gp", ym=.2, step=.2, cols=c(1,rgb(0,0,1,0.5)), zoom=NULL,
                            legend=TRUE, main="Density of low intensities"){
  ym0=min(y[y>0]) # ym=ym0 si fixe,  ym=param optimisé sinon
  tab=seq(ym0-step/2,max(y)+step/2,by=step)
  pemp=table(cut(y,tab))/length(y)
  tab2=round(seq(ym0,max(y),by=step),10)
  pth=dtcG(tab2,theta=theta, name=name,ym=ym, step=step)
  p=rbind(pemp,pth)
  colnames(p)=tab2
  if (all(zoom==range(y))){w=1:min(20,ncol(p))}else{w=which(tab2<=zoom[2] & tab2>=zoom[1])}
  barplot(p[,w],beside=TRUE, main=main, col=cols,ylab="p(Y=y)",xlab="y")
  if (legend) legend("topright",c("Empirical",paste("Best AIC tcG model:",name)),bty="n",text.col=cols)
}
