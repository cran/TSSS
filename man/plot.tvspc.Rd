\name{plot.tvspc}
\alias{plot.tvspc}

\title{Plot Evolutionary Power Spectra Obtained by Time Varying AR Model}
 
\description{
  Plot evolutionary power spectra obtained by time varying AR model returned
  by \code{\link{tvspc}}.
}
\usage{
\method{plot}{tvspc}(x, theta = 0, phi = 15, expand = 1, col = "lightblue",
        ticktype= "detail", \dots)
}
\arguments{
  \item{x}{an object of class \code{"tvspc"}.}
  \item{theta, phi, expand, col, ticktype}{graphical parameters in perspective
        plot \code{\link[graphics]{persp}}.}
  \item{\dots}{further graphical parameters may also be supplied as arguments.}
}
\examples{
# seismic data
data(MYE1F)
z <- tvar(MYE1F, trend.order = 2, ar.order = 8, span = 20,
          outlier = c(630, 1026), tau2.ini = 6.6e-06, delta = 1.0e-06)
spec <- tvspc(z$arcoef, z$sigma2)
plot(spec, theta = 30, phi = 40, expand = 0.5)
}
\keyword{ts}
