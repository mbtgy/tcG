#' Log-likelihood function of a meta-Gaussian model
#'
#'
#' @param discrete boolean. \code{TRUE} for the discrete likelihhod, \code{FALSE} for the continuous one
#' @param name name of the anamorphosis, one of \code{'power'},\code{'power-exp'}, \code{'quadratic-power'}, \code{'gp'}
#'
#' @importFrom stats pnorm
#' @return The log-likelihood function of meta-Gaussian model with anamorphosis \code{name}.
#' @export
#'
#' @examples
#' ll.choice(TRUE, "power")

ll.choice = function(discrete, name="power-exp"){
  if (discrete){
    if (name!="gp"){
      l=function(theta,y, name, ym=theta[length(theta)], step=0.2, sign=-1){
        p=sum(y==0)
        y = y[y>0]
        out=sign*(p*log(pnorm(-theta[1]))+
                    sum(log(pnorm(psi1(theta,y+step,ym,name)-theta[1])-
                              pnorm(psi1(theta,y,ym,name)-theta[1]))))
        if (is.na(out) | out==Inf) out=1e50
        return(out)
      }
    }else{
      l=function(theta,y, name, ym=theta[length(theta)], step=0.2, sign=-1){
        cst=pnorm(sqrt(abs(-1/(min(theta[3]*theta[4],0))))-theta[1])
        ysup=ym+theta[2]*(exp(-1)/max(c(-theta[3]*theta[4],0)))^(1/2/theta[3])
        p=sum(y==0)
        y = y[y>0]
        if(max(y)>ysup){out=1e50}else{
          out=sign*(p*log(pnorm(-theta[1])/cst)+
                      sum(log(pnorm(psi1(theta,pmin(y+step,ysup),ym,name)-theta[1])/cst-
                                pnorm(psi1(theta,pmax(y,ym),ym,name)-theta[1])/cst))) # pmax gère le cas ym paramètre
        }
        if (is.na(out) | out==Inf) out=1e50
        return(out)
      }
    }
  }else{
    l=function(theta,y, name, ym=theta[length(theta)], step=0.2, sign=-1){
      p=sum(y==0)
      y=y[y>0]
      return(sign*(p*log(pnorm(-theta[1]))+sum(log(dtcG(y,theta,name,ym,step)))))
    }
  }
  return(l)
}
