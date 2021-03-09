#' Fitting meta-Gaussian model(s)
#'
#' This function allows fitting 4 meta-Gaussian models.
#' The data should only contain positive values and have a discrete component in zero.
#'
#' @param y vector of data
#' @param name name of the anamorphosis that will be fitted, a subset of \code{c('gp', 'power', 'power-exp', 'quadratic-power')}
#' @param init initialisation parameters, a list named with the \code{name} argument
#' @param ym minimal value that can be observed
#' @param step discretization step
#' @param plots boolean, should plots be drawn?
#' @param bootstrap boolean, should boostrap replicates be fitted?
#' @param R number of bootstrap replicates
#' @param filename if you want to save the plots, a character for your file's name
#' @param ncpus number of cores
#' @param ym.param boolean, should \code{ym} be a parameter?
#'
#' @importFrom stats optim quantile
#' @importFrom grDevices rgb png dev.off
#' @importFrom graphics hist lines legend
#' @importFrom boot boot
#' @importFrom Rdpack reprompt
#'
#' @details A meta-Gaussian model is a variable \eqn{Y} defined as
#'
#' \eqn{Y = 0 * }\code{I}\eqn{_{X<0} + \psi(X) * }\code{I}\eqn{_{X\ge0}}
#'
#' with \eqn{X~N(\mu,1)} and where \code{I} is the indicator function, equals to 1 if the condition is true and 0 else.
#'
#' The available fun are the following:
#'   - \code{"gp"} for GP meta-Gaussian as in \insertCite{boutigny2021modelling;textual}{tcG}:
#'      \eqn{\psi(x) = ym+ \sigma x^(1/\alpha) exp(\xi x^2/2)}
#'
#'   - \code{"power"} as in \insertCite{bardossy1992space;textual}{tcG}: \eqn{\psi(x) = ym + \sigma x^(1/\alpha)}
#'
#'   - \code{"power-exp"} as in \insertCite{allard2015disaggregating;textual}{tcG}: \eqn{\psi(x) = ym + \sigma2 (exp(\sigma1 * x^(1/\alpha))-1)}
#'
#'   - \code{"quadratic-power"} as in \insertCite{allcroft2003latent;textual}{tcG}: \eqn{\psi(x) = ym + \sigma1 * x^(1/\alpha) + \sigma2 * x^(2/\alpha)}
#'
#' @return A list of four elements containing the results of the fitted model(s).
#' The first three elements are lists named by the fitted model(s) (argument \code{name}).
#' \code{par} gives the fitted parameters.
#' \code{AIC} gives the Akaike criterion.
#' \code{fit.boot} gives the parameters fitted with the bootstrap replicates.
#' Finally \code{for.plots} is a list that is used by the function \code{res.plot}.
#'
#' @references{
#' \insertAllCited{}
#' }
#'
#' @export
#'
#' @examples
#'\dontrun{
#' data(rain, package = "ismev")
#'
#' # The histogram is not very revelant due to the fact \code{rain} are daily observations
#' res=tcG.fit(rain, name=c("power", "gp"),
#'             init=list("power"=c(0,2,2), "gp"=c(0,2,0.5,0.5)),
#'             ym=0.3, step=0.1, R=50)
#'}

tcG.fit = function(y, name, init, ym=0, step=0, plots=TRUE, bootstrap=TRUE, R=500,
                   filename=NULL, ncpus=parallel::detectCores(), ym.param=FALSE){
  y=y[!is.na(y)]
  p=c(1:length(y))/(length(y) + 1)
  x=seq(0, max(y, na.rm = TRUE), by = 0.05)
  k=list("power"=3,"power-exp"=4,"quadratic-power"=4,"gp"=4)
  if (ym.param) lapply(k,function(x) x=x+1)

  AIC=par=fit.boot=q.L=q.U=qfit=dfit=list()
  for (nm in name){
    ll=ll.choice(step!=0,nm)
    print(nm)
    if (!ym.param){opt=optim(init[[nm]], fn=ll, y=y,ym=ym, step=step, name=nm)}else{
      opt=optim(init[[nm]], fn=ll, y=y, step=step, name=nm)
      ym=opt$par[k[[nm]]]}

    par[[nm]] = opt$par
    AIC[[nm]] = 2*opt$value + 2*k[[nm]]

    if (bootstrap){
      if (!ym.param){
        fit.boot[[nm]] = boot(data = y,function(data, original,ym, step, name, par){
          optim(par,fn=ll,y=data[original],ym=ym,step=step,name=name,control=list(maxit=2000))$par
        }, R=R, parallel="multicore",ncpus=ncpus,ym=ym, step=step, name=nm, par=par[[nm]])
      }else{
        fit.boot[[nm]] = boot(data = y,function(data, original, step, name, par){
          optim(par,fn=ll,y=data[original],step=step,name=name,control=list(maxit=2000))$par
        }, R=R, parallel="multicore",ncpus=ncpus, step=step, name=nm, par=par[[nm]])
      }

      q.boot = mapply(FUN=qtcG, p=list(c(1:length(y))/(length(y) + 1)),name=nm,ym=ym, step=step,
                      theta=split(t(fit.boot[[nm]]$t),rep(1:R, each=k[[nm]])))
      q.L[[nm]] = apply(q.boot,1,quantile,0.025,na.rm=TRUE)
      q.U[[nm]] = apply(q.boot,1,quantile,0.975,na.rm=TRUE)
      fit.boot[[nm]] = fit.boot[[nm]]$t
    }
    qfit[[nm]] = qtcG(p, name=nm, theta=par[[nm]],ym=ym, step=step)
    dfit[[nm]] = dtcG(x, name=nm, theta=par[[nm]],ym=ym, step=step)
  }

  if (plots){
    cols=cbind(c(rgb(0,0,1), rgb(1,0,0), rgb(0,1,0), rgb(1,0,1)),
               c(rgb(0,0,1,0.5), rgb(1,0,0,0.5), rgb(0,1,0,0.5), rgb(1,0,1,0.5)),
               c(rgb(0,0,1,0.1), rgb(1,0,0,0.1), rgb(0,1,0,0.1), rgb(1,0,1,0.1)))
    w=which.min(AIC)
    nm=name[w]
    if (ym.param) ym=par[[nm]][k[[nm]]]
    if (!is.null(filename)) png(filename, 1200,600)
    par(mfrow=c(1,2))
    if (step==0){
      hist(y[y>0], breaks=c(0:30)*max(y+.1)/30, freq=FALSE, xlim=range(x), main="Density of y>0",
           xlab="Precipitation [mm]", ylab="Density", col="lightgrey", cex.main=1)
      lines(x[x>=ym], dfit[[w]][x>=ym]/(1-pnorm(-par[[nm]][1])), col=cols[w,1])
      legend("topright", paste("Best AIC tcG model:",nm), text.col=cols[w,1], bty="n")
    }else{
      tcG.barplot(y,theta=par[[nm]],name=nm,ym=ym,step=step,cols=c(rgb(0,0,0),cols[w,2]), zoom=range(y,na.rm=T))
    }
    rg=range(y, qfit, q.L, q.U, na.rm=TRUE)
    tcG.qqplot(y, names=name, qfit=qfit, q.L=q.L, q.U=q.U, cols=cols, zoomx=rg, zoomy=rg)
    if (!is.null(filename)) dev.off()
  }
  return(list(par=par, fit.boot=fit.boot, AIC=AIC, for.plots=list(dfit=dfit,qfit=qfit,q.L=q.L,q.U=q.U)))
}
