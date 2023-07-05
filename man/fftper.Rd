\name{fftper}
\alias{fftper}

\title{Compute a Periodogram via FFT}
\description{
  Compute a periodogram of the univariate time series via FFT.}
\usage{
fftper(y, window = 1, plot = TRUE, \dots)
}
\arguments{
  \item{y}{a univariate time series.}
  \item{window}{smoothing window type. (0: box-car, 1: Hanning, 2: Hamming)}
  \item{plot}{logical. If \code{TRUE} (default), smoothed (log-)periodogram is
              plotted.}
  \item{\dots}{graphical arguments passed to \code{\link{plot.spg}}.}
}
\details{
  \tabular{rll}{
    Hanning Window : \tab \eqn{W_0} = 0.5 \tab \eqn{W_1} = 0.25 \cr
    Hamming Window : \tab \eqn{W_0} = 0.54 \tab \eqn{W_1} = 0.23}
}
\value{
  An object of class \code{"spg"}, which is a list with the following components:
  \item{period}{periodogram.}
  \item{smoothed.period}{smoothed periodogram. If there is not a negative
    number, logarithm of smoothed periodogram.}
  \item{log.scale}{logical. If \code{TRUE} \code{smoothed.period} is logarithm
    of smoothed periodogram.}
  \item{tsname}{the name of the univariate time series \code{y}.}
}

\note{
  We assume that the length \eqn{N} of the input time series \code{y} is a power
  of 2. If \eqn{N} is not a power of 2, calculate using the FFT by appending 0's
  behind the data \code{y}. 
}
\references{
  Kitagawa, G. (2020)
  \emph{Introduction to Time Series Modeling with Applications in R}.
  Chapman & Hall/CRC.
}
\examples{
# Yaw rate, rolling, pitching and rudder angle of a ship
data(HAKUSAN)
YawRate <- HAKUSAN[, 1]
fftper(YawRate, window = 0)
}

\keyword{ts}

