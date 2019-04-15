\name{lsar.chgpt}
\alias{lsar.chgpt}

\title{Estimation of the Change Point}
\description{
  Precisely estimate a change point of subinterval for locally stationary AR
  model.
}
\usage{
lsar.chgpt(y, max.arorder = 20, subinterval, candidate, plot = TRUE, \dots)
}
\arguments{
  \item{y}{a univariate time series.}
  \item{max.arorder}{highest order of AR model.}
  \item{subinterval}{a vector of the form \code{c(n0, ne)} which gives a start
    and end point of time interval used for model fitting.}
  \item{candidate}{a vector of the form \code{c(n1, n2)} which gives minimum
    and maximum for change point.

    \code{n0+2k} < \code{n1} < \code{n2+k} < \code{ne},

    (\code{k} is \code{max.arorder})}
  \item{plot}{logical. If \code{TRUE} (default), \code{y[n0:ne]} and
    '\code{aic}' are plotted.}
  \item{\dots}{further arguments to be passed to \code{plot.chgpt}.}
}
\value{
  An object of class \code{"chgpt"}, which is a list with the following
  elements:
  \item{aic}{AICs of the AR model fitted on \code{[n1, n2]}.}
  \item{aicmin}{minimum AIC.}
  \item{change.point}{a change point.}
  \item{subint}{original sub-interval data and information.} 
}
\references{
  Kitagawa, G. (2010)
  \emph{Introduction to Time Series Modeling}. Chapman & Hall/CRC.
}
\examples{
# seismic data
data(MYE1F)
lsar.chgpt(MYE1F, max.arorder = 10, subinterval = c(200, 1000),
           candidate = c(400, 800))

lsar.chgpt(MYE1F, max.arorder = 10, subinterval = c(600, 1400),
           candidate = c(800, 1200))
}

\keyword{ts}
