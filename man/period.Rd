\name{period}
\alias{period}

\title{Compute a Periodogram}
\description{
  Compute a periodogram of the univariate time series.}
\usage{
period(y, window = 1, minmax = c(-1.0e+30, 1.0e+30), plot = TRUE, \dots)
}
\arguments{
  \item{y}{a univariate time series.}
  \item{window}{smoothing window type. (0: box-car, 1: Hanning, 2: Hamming)}
  \item{minmax}{bound for outliers in low side and high side.}
  \item{plot}{logical. If \code{TRUE} (default), smoothed periodogram is plotted.}
  \item{\dots}{further arguments to be passed to \code{\link{plot.spg}}.}
}
\details{
  \tabular{rll}{
    Hanning Window : \tab \eqn{W_0} = 0.5 \tab \eqn{W_1} = 0.25 \cr
    Hamming Window : \tab \eqn{W_0} = 0.54 \tab \eqn{W_1} = 0.23}
}
\value{
  An object of class \code{"spg"}, which is a list with the following
  elements:
  \item{period}{periodogram (raw spectrum).}
  \item{smoothed.period}{smoothed periodogram. If there is not a negative number,
    logarithm of smoothed periodogram.}
  \item{log.scale}{if \code{TRUE} "smooth the periodogram on log scale.}
  \item{tsname}{the name of the univariate time series \code{y}.}
}

\references{
  Kitagawa, G. (2010)
  \emph{Introduction to Time Series Modeling}. Chapman & Hall/CRC.
}
\examples{
# BLSALLFOOD data
data(BLSALLFOOD)
period(BLSALLFOOD)

# seismic Data
data(MYE1F)
period(MYE1F)
}

\keyword{ts}
