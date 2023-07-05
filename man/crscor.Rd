\name{crscor}
\alias{crscor}

\title{Cross-Covariance and Cross-Correlation}
\description{
  Compute cross-covariance and cross-correlation functions of the multivariate
  time series.
}
\usage{
crscor(y, lag = NULL, outmin = NULL, outmax = NULL, plot = TRUE, \dots)
}
\arguments{
  \item{y}{a multivariate time series.}
  \item{lag}{maximum lag. Default is \eqn{2 \sqrt{n}}{2*sqrt(n)}, where \eqn{n}
    is the length of the time series \code{y}.}
  \item{outmin}{bound for outliers in low side. A default value is -1.0e+30 for
    each dimension.}
  \item{outmax}{bound for outliers in high side. A default value is 1.0e+30 for
    each dimension.}
  \item{plot}{logical. If \code{TRUE} (default), cross-correlations are plotted.}
  \item{\dots}{graphical arguments passed to the \code{plot} method.}
}
\value{
  An object of class \code{"crscor"} which has a \code{plot} method. This is a
  list with the following components:
  \item{cov}{cross-covariances.}
  \item{cor}{cross-correlations.}
  \item{mean}{mean vector.}
}
\references{
  Kitagawa, G. (2020)
  \emph{Introduction to Time Series Modeling with Applications in R}.
  Chapman & Hall/CRC.
}
\examples{
# Yaw rate, rolling, pitching and rudder angle of a ship
data(HAKUSAN)
y <- as.matrix(HAKUSAN[, 2:4])   # Rolling, Pitching, Rudder
crscor(y, lag = 50)

# The groundwater level and the atmospheric pressure
data(Haibara)
crscor(Haibara, lag = 50)
}

\keyword{ts}
