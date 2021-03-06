\name{unicor}
\alias{unicor}

\title{Autocovariance and Autocorrelation}
\description{
  Compute autocovariance and autocorrelation function of the univariate time
 series.}
\usage{
unicor(y, lag = NULL, minmax = c(-1.0e+30, 1.0e+30), plot = TRUE, \dots)
}
\arguments{
  \item{y}{a univariate time series.}
  \item{lag}{maximum lag. Default is \eqn{2 \sqrt{n}}{2*sqrt(n)}, where \eqn{n}
    is the length of the time series \code{y}.}
  \item{minmax}{thresholds for outliers in low side and high side.}
  \item{plot}{logical. If \code{TRUE} (default), autocorrelations are plotted.}
  \item{\dots}{graphical arguments passed to the \code{plot} method.}
}
\value{
  An object of class \code{"unicor"} which has a \code{plot} method. This is a
  list with the following components:
  \item{acov}{autocovariances.}
  \item{acor}{autocorrelations.}
  \item{acov.err}{error bound for autocovariances.}
  \item{acor.err}{error bound for autocorrelations.}
  \item{mean}{mean of \code{y}.}
  \item{tsname}{the name of the univariate time series \code{y}.}
}
\references{
  Kitagawa, G. (2020)
  \emph{Introduction to Time Series Modeling with Applications in R}.
  Chapman & Hall/CRC.
}
\examples{
# Yaw rate, rolling, pitching and rudder angle of a ship
data(HAKUSAN)
Yawrate <- HAKUSAN[, 1]
unicor(Yawrate, lag = 50)

# seismic data
data(MYE1F)
unicor(MYE1F, lag = 50)
}

\keyword{ts}
