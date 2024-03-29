% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_tcG.R
\name{tcG.fit}
\alias{tcG.fit}
\title{Fitting meta-Gaussian model(s)}
\usage{
tcG.fit(
  y,
  name,
  init,
  ym = 0,
  step = 0,
  plots = TRUE,
  bootstrap = TRUE,
  R = 500,
  filename = NULL,
  ncpus = parallel::detectCores(),
  ym.param = FALSE
)
}
\arguments{
\item{y}{vector of data}

\item{name}{name of the anamorphosis that will be fitted, a subset of \code{c('gp', 'power', 'power-exp', 'quadratic-power')}}

\item{init}{initialisation parameters, a list named with the \code{name} argument}

\item{ym}{minimal value that can be observed}

\item{step}{discretization step}

\item{plots}{boolean, should plots be drawn?}

\item{bootstrap}{boolean, should boostrap replicates be fitted?}

\item{R}{number of bootstrap replicates}

\item{filename}{if you want to save the plots, a character for your file's name}

\item{ncpus}{number of cores}

\item{ym.param}{boolean, should \code{ym} be a parameter?}
}
\value{
A list of four elements containing the results of the fitted model(s).
The first three elements are lists named by the fitted model(s) (argument \code{name}).
\code{par} gives the fitted parameters.
\code{AIC} gives the Akaike criterion.
\code{fit.boot} gives the parameters fitted with the bootstrap replicates.
Finally \code{for.plots} is a list that is used by the function \code{res.plot}.
}
\description{
This function allows fitting 4 meta-Gaussian models.
The data should only contain positive values and have a discrete component in zero.
}
\details{
A meta-Gaussian model is a variable \eqn{Y} defined as

\eqn{Y = 0 * }\code{I}\eqn{_{X<0} + \psi(X) * }\code{I}\eqn{_{X\ge0}}

with \eqn{X~N(\mu,1)} and where \code{I} is the indicator function, equals to 1 if the condition is true and 0 else.

The available fun are the following:
\itemize{
\item \code{"gp"} for GP meta-Gaussian as in \insertCite{boutigny2021modelling;textual}{tcG}:
\eqn{\psi(x) = ym+ \sigma x^(1/\alpha) exp(\xi x^2/2)}
\item \code{"power"} as in \insertCite{bardossy1992space;textual}{tcG}: \eqn{\psi(x) = ym + \sigma x^(1/\alpha)}
\item \code{"power-exp"} as in \insertCite{allard2015disaggregating;textual}{tcG}: \eqn{\psi(x) = ym + \sigma1 (exp(\sigma2 * x^(1/\alpha))-1)}
\item \code{"quadratic-power"} as in \insertCite{allcroft2003latent;textual}{tcG}: \eqn{\psi(x) = ym + \sigma1 * x^(1/\alpha) + \sigma2 * x^(2/\alpha)}
}
}
\examples{
\dontrun{
data(rain, package = "ismev")

# The histogram is not very revelant due to the fact \code{rain} are daily observations
res=tcG.fit(rain, name=c("power", "gp"),
            init=list("power"=c(0,2,2), "gp"=c(0,2,0.5,0.5)),
            ym=0.3, step=0.1, R=50)
}
}
\references{
{
\insertAllCited{}
}
}
