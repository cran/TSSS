\name{plot.spg}
\alias{plot.spg}
\title{Plot Smoothed Periodogram}
\description{
  Plot smoothed periodogram or logarithm of smoothed periodogram.
}
\usage{
\method{plot}{spg}(x, type = "vl", \dots)
}
\arguments{
  \item{x}{an object of class \code{"spg"} as returned by \code{\link{period}} 
    and \code{\link{fftper}}.}
  \item{type}{type of plot. ("l": lines, "vl" : vertical lines)}
  \item{\dots}{further graphical parameters may also be supplied as arguments.}
}
\keyword{ts}
