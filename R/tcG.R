#' Meta-Gaussian models distributions
#'
#' Density, distribution function, quantile function and random generation of a meta-Gaussian distribution with parameters theta.
#'
#'
#' @param y vector of quantiles (rainfall)
#' @param theta model parameters
#' @param name name of the anamorphosis, one of \code{'power'},\code{'power-exp'}, \code{'quadratic-power'}, \code{'gp'}
#' @param ym minimal value that can be observed
#' @param step discretization step
#'
#' @importFrom VGAM lambertW
#' @importFrom stats pnorm dnorm qnorm rnorm
#'
#' @details Four meta-Gaussian models are available (name argument): \code{'power'},\code{'power-exp'}, \code{'quadratic-power'}, \code{'gp'}
#'
#'
#' @return \code{dtcG} gives the density, \code{ptcG} gives the distribution function, \code{qtcG} gives the quantile function and \code{rtcG} generates random deviates.
#'
#' @examples
#' theta=c(-1,0.2,0.9,0.1)
#' s=seq(0,2,0.01)
#'
#' sim = rtcG(1e4, theta, "gp")
#' d = dtcG(s, theta, "gp")
#' hist(sim, probability=TRUE)
#' lines(s, d, col=2)



#' @export
dtcG = function(y,theta, name, ym=0, step=0){
  if (name=="gp"){
    cst=pnorm(sqrt(abs(1/min(theta[3]*theta[4],0)))-theta[1]) #P(X<xsup), cste normalisation pour troncation
    ysup=ym+theta[2]*(exp(-1)/max(-theta[3]*theta[4],0))^(1/2/theta[3]) # borne sup du domaine
  }else{cst=1;ysup=Inf}
  d=y
  nuls=(y==0 | y<ym | y>=ysup)
  d[nuls]=0 #p(y=y)=0 si y<ym et y>=ysup, except y=0
  d[y==0]=pnorm(-theta[1])/cst # cas y=0

  if (step==0){
    y=y[!nuls]
    dpsi1=switch(name,
                 "power"=(y-ym)^(1/theta[3]-1)/(theta[3]*theta[2]^(1/theta[3])),
                 "power-exp"=(log((y-ym)/theta[3]+1))^(1/theta[4]-1)/(theta[4]*theta[2]^(1/theta[4])*(y-ym+theta[3])),
                 "quadratic-power"=(2*theta[3]/(sqrt(theta[2]^2+4*theta[3]*(y-ym))-theta[2]))^(1-1/theta[4])/(theta[4]*sqrt(theta[2]^2+4*theta[3]*(y-ym))),
                 "gp"=sqrt(theta[3]*lambertW(theta[3]*theta[4]*((y-ym)/theta[2])^(2*theta[3]))/theta[4])/((y-ym)*(1+lambertW(theta[3]*theta[4]*((y-ym)/theta[2])^(2*theta[3])))))
    d[!nuls]=dnorm(psi1(theta, y, ym, name),theta[1])*dpsi1/cst
  }else{
    d[!nuls]=(pnorm(psi1(theta,y[!nuls]+step,ym,name)-theta[1])-
                pnorm(psi1(theta,y[!nuls],ym,name)-theta[1]))/cst
    d[y!=round(step*round(y/step),10)]=0
  }
  return(d)
}

#' @describeIn dtcG probability function
#' @export
ptcG = function(y, theta, name, ym=0, step=0){
  p=y
  nuls=(y==0)
  if (name=="gp"){
    cst=pnorm(sqrt(abs(1/min(theta[3]*theta[4],0)))-theta[1])  #P(X<xsup), cste normalisation pour troncation
  }else{cst=1}
  p[nuls]=pnorm(0, theta[1])/cst
  p[!nuls] = pnorm(psi1(theta, p[!nuls]+step, ym, name), mean=theta[1])/cst
  if (name=="gp"){p[is.na(p)]=1} # cas où y>sup
  return(p)
}

#' @param p vector of probabilities
#' @describeIn dtcG quantile function
#' @export
qtcG = function(p, theta, name, ym=0, step=0){
  if(name=="gp"){
    cst=pnorm(sqrt(abs(1/min(theta[3]*theta[4],0)))-theta[1])  #P(X<xsup), cste normalisation pour troncation
    q=psi(theta, qnorm(p*cst, theta[1]), ym=ym, name=name)#-step
    q[is.na(q)]=ym+theta[2]*(exp(-1)/max(-theta[3]*theta[4],0))^(1/2/theta[3]) # cas p=1, xi<0
  }else{
    q=psi(theta, qnorm(p, theta[1]), ym=ym, name=name)#-step
  }
  if (step!=0) q=step*floor(q/step)
  q[p<pnorm(-theta[1])]=0
  return(q)
}

#' @param n number of observations to generate
#' @describeIn dtcG random generation function
#' @export
rtcG = function(n, theta, name, ym=0, step=0){
  y=x=rnorm(n, theta[1])
  y[x<=0]=0
  y[x>0]=psi(theta, y[y>0], ym, name)
  y=y[!is.na(y)] # cas xi<0
  if (step!=0) y=step*floor(y/step)
  return(y)
}
