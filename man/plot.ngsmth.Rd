\name{plot.ngsmth}
\alias{plot.ngsmth}
\title{Plot Smoothed Density Function}
\description{
  Plot the smoothed density function returned by \code{\link{ngsmth}}.
}
\usage{
\method{plot}{ngsmth}(x, type = c("trend", "smt"), theta = 0, phi = 15,
                      expand = 1, col = "lightblue", ticktype= "detail", \dots)
}
\arguments{
  \item{x}{an object of class \code{"ngsmth"}.}
  \item{type}{plotted values, either or both of "trend" and "smt".}
  \item{theta, phi, expand, col, ticktype}{graphical parameters in perspective
        plot \code{\link[graphics]{persp}}.}
  \item{\dots}{further graphical parameters may also be supplied as arguments.}
}
\keyword{ts}
