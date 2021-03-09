#' Plotting the results
#'
#' This function offers choices for plotting results of \code{fit.tcG} and \code{fit.extGP}.
#' Two types of plots are available: a density plot (histogram or barplot) and a QQ plot.
#' In particular the function allows plotting the results of both fit functions together on the QQ plot.
#'
#' @param res.tcG,res.extGP a list resulting respectively from \code{fit.tcG} and \code{fit.extGP}, at least one of them must be non-null
#' @param y vector of data (if both \code{res.tcG} and \code{res.extGP} are given, \code{y} must include the zeros)
#' @param ym minimal value that can be observed
#' @param step discretization step
#' @param zoom vector of length 2, range to zoom on
#' @param select for \code{res.tcG}, a the model(s) that should be plotted. A subset of \code{c('gp', 'power', 'power-exp', 'quadratic-power')}
#' @param choice A subset of \code{c('qqplot', 'dens')}. Which plots should be drawn?
#' @param legend boolean, should there be a legend?
#' @param ... graphical parameters such as \code{main}, \code{cex}, etc.
#'
#' @importFrom grDevices rgb
#' @importFrom graphics lines hist legend par
#' @importFrom stats pnorm
#' @importFrom methods hasArg
#'
#' @export
#'
#' @examples
#'\dontrun{
#'
#' # Comparing meta-Gaussian models
#'
#' data(rain, package = "ismev")
#'
#' res=tcG.fit(rain, name=c("power", "gp"),
#'             init=list("power"=c(0,2,2), "gp"=c(0,2,0.5,0.5)),
#'             ym=0.3, step=0.1, R=50, plots=FALSE)
#' res.plot(res.tcG=res, y=rain, ym=.3, step=.1, zoom=c(0,30), choice="qqplot")
#'
#' # Comparing the GP meta-Gaussian and the extended GP
#'
#' res2=extGP.fit(rain[rain>0], c(2,.5,.5), ym=.3, step=.1, R=50, plots=FALSE)
#' res.plot(res.tcG=res, res.extGP=res2, y=rain, ym=.3, step=.1, select="gp", choice="qqplot")
#'
#'}
res.plot = function(res.tcG=NULL, res.extGP=NULL, y, ym=0, step=0, zoom=range(y,na.rm=TRUE),
                    select=NULL, choice=c("dens","qqplot"), legend=TRUE,...){
  y=y[!is.na(y)]
  if (!is.null(select)) names=select else names=names(res.tcG$par)
  p=c(1:length(y))/(length(y) + 1)
  x=seq(ym, max(y, na.rm = TRUE), by = 0.05)
  cols = cbind(c(rgb(0,0,1), rgb(1,0,0), rgb(0,1,0), rgb(1,0,1), rgb(.6,.3,0)),
               c(rgb(0,0,1,0.5), rgb(1,0,0,0.5), rgb(0,1,0,0.5), rgb(1,0,1,0.5), rgb(.6,.3,0,.5)),
               c(rgb(0,0,1,0.1), rgb(1,0,0,0.1), rgb(0,1,0,0.1), rgb(1,0,1,0.1), rgb(.6,.3,0,.1)))
  par(mfrow=c(1,length(choice)))
  if ("dens" %in% choice){
    if (!is.null(res.tcG)){w=which.min(res.tcG$AIC[names]);nm=names[w]}
    if (step==0){ # histogram
      if (any(zoom!=range(y,na.rm=T))){
        m=max(zoom[1], min(y[y>0]))
        y2=y[y>=m & y<=zoom[2]]
        if (round(m,2)==0){p="]"}else{p='['}
        if (!is.null(res.tcG)){
          cst=ptcG(zoom[2],res.tcG$par[[nm]],nm,ym,step)-ptcG(m,res.tcG$par[[nm]],nm,ym,step)
        }else{cst=pextGP(zoom[2],res.extGP$par,ym)-pextGP(m,res.extGP$par,ym)}
        if (!hasArg(main)) main=bquote("Density of y" %in% .(paste(p,round(m,2),",",round(zoom[2],2),"]",sep="")))
      }else{
        y2=y
        if (!is.null(res.tcG)){cst=1-pnorm(-res.tcG$par[[nm]][1])}else{cst=1}
        if (!hasArg(main)) main="Density of y>0"
      }

      if (!hasArg(main)){hist(y2[y2>0], breaks=c(0:30)*max(y2+.1)/30, freq=FALSE, xlim=zoom,
                              xlab="Precipitation [mm]", ylab="Density", col="lightgrey",main=main,cex.main=1)
      }else{hist(y2[y2>0], breaks=c(0:30)*max(y2+.1)/30, freq=FALSE, xlim=zoom,
                 xlab="Precipitation [mm]", ylab="Density", col="lightgrey",cex.main=1,...)}

      if (!is.null(res.tcG)){lines(x, res.tcG$for.plots$dfit[[w]]/cst, col=cols[w,1])
      }else{lines(x, res.extGP$for.plots$dfit/cst, col=cols[5,1])}
      if (legend) legend("topright", paste("Best AIC tcG model:",nm), text.col=cols[w,1], bty="n")
    }else{ # barplot
      if (!is.null(res.tcG)){tcG.barplot(y,theta=res.tcG$par[[nm]],name=nm,ym=ym,step=step,zoom=zoom,cols=c(rgb(0,0,0),cols[w,2]), legend=legend, ...)
      }else{extGP.barplot(y[y>0],theta=res.extGP$par,ym=ym,step=step,c(1,cols[5,2]),zoom=zoom, legend=legend,...)}
    }
  }
  if ("qqplot" %in% choice){
    if (all(zoom == range(y))){
      zoomy=range(y, res.tcG$for.plots$qfit, res.tcG$for.plots$q.L, res.tcG$for.plots$q.U,
                  res.extGP$for.plots$qfit, res.extGP$for.plots$q.L, res.extGP$for.plots$q.U, na.rm=TRUE)
      zoomx=range(y)
    }else{zoomx=zoomy=zoom}

    if (!is.null(res.tcG)){
      tcG.qqplot(y,names=names,qfit=res.tcG$for.plots$qfit,q.L=res.tcG$for.plots$q.L,
                 q.U=res.tcG$for.plots$q.U,cols=cols,zoomx=zoomx,zoomy=zoomy, legend=legend,...)
    }
    if (!is.null(res.extGP)){
      extGP.qqplot(y[y>0],res.extGP$for.plots$qfit, res.extGP$for.plots$q.L,
                   res.extGP$for.plots$q.U, add=(!is.null(res.tcG)), zoomx=zoomx,zoomy=zoomy, legend=legend,...)
    }
  }
}
