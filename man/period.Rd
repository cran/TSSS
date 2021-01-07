\name{period}
\alias{period}

\title{Compute a Periodogram}
\description{
  Compute a periodogram of the univariate time series.}
\usage{
period(y, window = 1, lag = NULL, minmax = c(-1.0e+30, 1.0e+30),
       plot = TRUE, \dots)
}
\arguments{
  \item{y}{a univariate time series.}
  \item{window}{smoothing window type. (0: box-car, 1: Hanning, 2: Hamming)}
  \item{lag}{maximum lag of autocovariance. If \code{NULL} (default), \cr
    \code{window} = 0 : \code{lag} = \eqn{n} - 1, \cr
    \code{window} > 0 : \code{lag} = 2 \eqn{\sqrt{n}}{sqrt(n)}, \cr
    where \eqn{n} is the length of data.}
  \item{minmax}{bound for outliers in low side and high side.}
  \item{plot}{logical. If \code{TRUE} (default), smoothed periodogram is plotted.}
  \item{\dots}{graphical arguments passed to \code{\link{plot.spg}}.}
}
\details{
  \tabular{rll}{
    Hanning Window : \tab \eqn{W_0} = 0.5 \tab \eqn{W_1} = 0.25 \cr
    Hamming Window : \tab \eqn{W_0} = 0.54 \tab \eqn{W_1} = 0.23}
}
\value{
  An object of class \code{"spg"}, which is a list with the following components:
  \item{period}{periodogram(or raw spectrum).}
  \item{smoothed.period}{smoothed log-periodogram. Smoothed periodogram is given
    if there is a negative value in the smoothed periodogram.}
  \item{log.scale}{if \code{TRUE} "smooth the periodogram on log scale.}
  \item{tsname}{the name of the univariate time series \code{y}.}
}
\references{
  Kitagawa, G. (2020)
  \emph{Introduction to Time Series Modeling with Applications in R}.
  Chapman & Hall/CRC.
}
\examples{
## BLSALLFOOD data
data(BLSALLFOOD)
period(BLSALLFOOD)

## seismic Data
data(MYE1F)

# smoothed periodogram 
period(MYE1F)

# periodogram
period(MYE1F, window = 0)

# raw spectrum
period(MYE1F, window = 0, lag = 200)

# Hamming window
period(MYE1F, window = 2)
}

\keyword{ts}
