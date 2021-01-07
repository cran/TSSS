\name{plot.tvspc}
\alias{plot.tvspc}

\title{Plot Evolutionary Power Spectra Obtained by Time Varying AR Model}
 
\description{
  Plot evolutionary power spectra obtained by time varying AR model returned
  by \code{\link{tvspc}}.
}
\usage{
\method{plot}{tvspc}(x, tvv = NULL, dx = 2, dy = 0.25, \dots)
}
\arguments{
  \item{x}{an object of class \code{"tvspc"}.}
  \item{tvv}{time varying variance as returned by \code{\link{tvvar}}.}
  \item{dx}{step width for the X axis.}
  \item{dy}{step width for the Y axis.}
  \item{\dots}{further graphical parameters may also be supplied as arguments.}
}
\examples{
# seismic data
data(MYE1F)
v <- tvvar(MYE1F, trend.order = 2, tau2.ini = 6.6e-06, delta = 1.0e-06,
           plot = FALSE )

z <- tvar(v$nordata, trend.order = 2, ar.order = 8, span = 20,
          outlier = c(630, 1026), tau2.ini = 6.6e-06, delta = 1.0e-06,
          plot = FALSE)

spec <- tvspc(z$arcoef, z$sigma2, span = 20, nf = 400)
plot(spec, tvv = v$tvv, dx = 2, dy = 0.10)
}

\keyword{ts}
