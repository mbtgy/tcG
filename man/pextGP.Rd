% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extGP.R
\name{pextGP}
\alias{pextGP}
\alias{qextGP}
\alias{dextGP}
\title{Extended Generalized Pareto distribution}
\usage{
pextGP(y, theta, ym = 0)

qextGP(p, theta, ym = 0, step = 0)

dextGP(y, theta, ym = 0, step = 0)
}
\arguments{
\item{y}{vector of quantiles (rainfall)}

\item{theta}{model parameters (sigma, alpha, xi)}

\item{ym}{minimal value that can be observed}

\item{p}{vector of probabilities}

\item{step}{discretization step}
}
\value{
\code{dextGP} gives the density, \code{pextGP} gives the distribution function, \code{qextGP} gives the quantile function and \code{rextGP} generates random deviates.
}
\description{
Density, distribution function, quantile function and random generation of an extended Generalised Pareto distribution with parameters theta.
}
\section{Functions}{
\itemize{
\item \code{qextGP}: quantile function

\item \code{dextGP}: density function
}}

