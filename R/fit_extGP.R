#' Fitting the extended GP model
#'
#' @param y vector of positive data
#' @param init initialisation parameters, a vector of length 3
#' @param ym minimal value that can be observed
#' @param step discretization step
#' @param plots boolean, should plots be drawn?
#' @param bootstrap boolean, should boostrap replicates be fitted?
#' @param R number of bootstrap replicates
#' @param filename if you want to save the plots, a character for your file's name
#' @param ncpus number of cores
#' @param ym.param boolean, should \code{ym} be a parameter?
#'
#' @details The extended GP model \insertCite{naveau2016modeling}{tcG} is define for positive rainfall as
#'
#' \eqn{Y+ = ym + \sigma H_\xi^{-1}(U^{1/\alpha})},
#'
#' with \eqn{H_\xi} the cdf of a GPD.
#'
#' @return A list of four elements containing the results of the fitted model.
#' \code{par} gives the fitted parameters.
#' \code{AIC} gives the Akaike criterion.
#' \code{fit.boot} gives the parameters fitted with the bootstrap replicates.
#' Finally \code{for.plots} is a list that is used by the function \code{res.plot}.
#' @export
#'
#' @references{
#' \insertAllCited{}
#' }
#'
#'
#' @examples
#'\dontrun{
#' res2=extGP.fit(rain[rain>0], c(2,.5,.5), ym=.3, step=.1, R=50)
#'}


extGP.fit = function(y, init, ym=0, step=0, bootstrap=TRUE, R=500,plots=TRUE,
                     filename=NULL, ncpus=parallel::detectCores(), ym.param=FALSE){
  y=y[which(y>0)]
  p=c(1:length(y))/(length(y) + 1)
  x=seq(0, max(y, na.rm = TRUE), by = 0.05)
  k=3

  ll=function(theta,y,ym,step){-sum(log(pextGP(y+step,theta,ym)-pextGP(pmax(y,ym),theta,ym)))}
  opt=optim(init, fn=ll, y=y,ym=ym, step=step)
  par = opt$par
  AIC = 2*opt$value + 2*k

  fit.boot=q.L=q.U=NULL
  if (bootstrap){
    fit.boot = boot::boot(data = y,function(data, original,ym, step, par){
      optim(par,fn=ll,y=data[original],ym=ym,step=step,control=list(maxit=2000))$par
    }, R=R, parallel="multicore",ncpus=ncpus,ym=ym, step=step, par=par)
    q.boot = mapply(FUN=qextGP, p=list(c(1:length(y))/(length(y) + 1)),ym=ym, step=step,
                    theta=split(t(fit.boot$t),rep(1:R, each=k)))
    q.L = apply(q.boot,1,quantile,0.025,na.rm=TRUE)
    q.U = apply(q.boot,1,quantile,0.975,na.rm=TRUE)
    fit.boot = fit.boot$t
  }
  qfit = qextGP(p, theta=par,ym=ym, step=step)
  dfit = dextGP(x, theta=par,ym=ym, step=step)


  if (plots){
    cols=c(rgb(.6,.3,0), rgb(.6,.3,0,0.5),rgb(.6,.3,0,0.1))
    if (!is.null(filename)) png(filename, 1200,600)
    par(mfrow=c(1,2))
    if (step==0){
      hist(y, breaks=c(0:30)*max(y+.1)/30, freq=FALSE, xlim=range(x), main="Density of y>0",
           xlab="Precipitation [mm]", ylab="Density", col="lightgrey")
      lines(x[x>=ym], dfit[x>=ym], col=cols[1])
      legend("topright", "extGP", text.col=cols[1], bty="n")
    }else{
      extGP.barplot(y[y>0],theta=par,ym=ym,step=step,cols=c(rgb(0,0,0),cols[2]), zoom=range(y,na.rm=T))
    }
    rg=range(y, qfit, q.L, q.U, na.rm=TRUE)
    extGP.qqplot(y[y>0], qfit=qfit, q.L=q.L, q.U=q.U, cols=cols, zoomx=rg, zoomy=rg)
    if (!is.null(filename)) dev.off()
  }
  return(list(par=par, fit.boot=fit.boot, AIC=AIC, for.plots=list(dfit=dfit,qfit=qfit,q.L=q.L,q.U=q.U)))
}
