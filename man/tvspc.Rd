\name{tvspc}
\alias{tvspc}

\title{Evolutionary Power Spectra by Time Varying AR Model}
\description{
  Estimate evolutionary power spectra by time varying AR model.
}
\usage{
tvspc(arcoef, sigma2, var = NULL, span = 20, nf = 200)
}
\arguments{
  \item{arcoef}{time varying AR coefficients.}
  \item{sigma2}{variance of the observational noise.}
  \item{var}{time varying variance.}
  \item{span}{local stationary span.}
  \item{nf}{number of frequencies in evaluating spectrum.}
}
\value{
  return an object of class \code{"tvspc"}.
}
\references{
  Kitagawa, G. (2010)
  \emph{Introduction to Time Series Modeling}. Chapman & Hall/CRC.

  Kitagawa, G. and Gersch, W. (1996)
  \emph{Smoothness Priors Analysis of Time Series}. Lecture Notes in Statistics,
  No.116, Springer-Verlag.

  Kitagawa, G. and Gersch, W. (1985)
  \emph{A smoothness priors time varying AR coefficient modeling 
  of nonstationary time series}. IEEE trans. on Automatic Control, AC-30, 48-56.
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
