#' Barplot for the extended GP model
#'
#'
#' @param y vector of data
#' @param theta vector of model parameters
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
extGP.barplot = function(y, theta, ym=.2, step=.2, cols=c(1,rgb(.6,.3,0,0.5)),
                         zoom=NULL, legend=TRUE, main="Density of positive low intensities"){
  y=y[y>0]
  ym0=min(y) # ym=ym0 si fixe,  ym=param optimisé sinon
  tab=seq(ym0-step/2,max(y)+step/2,by=step)
  pemp=table(cut(y,tab))/length(y)
  tab2=round(seq(ym0,max(y),by=step),10)
  pth=dextGP(tab2,theta=theta,ym=ym, step=step)
  p=rbind(pemp,pth)
  colnames(p)=tab2
  if (all(zoom==range(y))){w=1:min(ncol(p),20)}else{w=which(tab2<=zoom[2] & tab2>=zoom[1])}
  barplot(p[,w],beside=TRUE,main=main,col=cols,ylab="p(Y=y)",xlab="y",cex.main=1)
  if (legend) legend("topright",c("Empirical","Extended GP"),cex=.8,bty="n",text.col=cols)
}
