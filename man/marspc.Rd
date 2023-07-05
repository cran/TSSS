\name{marspc}
\alias{marspc}

\title{Cross Spectra and Power Contribution}
\description{
  Compute cross spectra, coherency and power contribution.
}
\usage{
marspc(arcoef, v, plot = TRUE, \dots)
}
\arguments{
  \item{arcoef}{AR coefficient matrices.}
  \item{v}{innovation variance matrix.}
  \item{plot}{logical. If \code{TRUE} (default), cross spectra, coherency and
    power contribution are plotted.}
  \item{\dots}{graphical arguments passed to the \code{plot} method.}
}
\value{
  An object of class \code{"marspc"} which has a \code{plot} method. This is a
  list with the following components:
  \item{spec}{cross spectra.}
  \item{amp}{amplitude spectra.}
  \item{phase}{phase spectra.}
  \item{coh}{simple coherency.}
  \item{power}{decomposition of power spectra.}
  \item{rpower}{relative power contribution.}
}
\references{
  Kitagawa, G. (2020)
  \emph{Introduction to Time Series Modeling with Applications in R}.
  Chapman & Hall/CRC.
}
\examples{
# Yaw rate, rolling, pitching and rudder angle of a ship
data(HAKUSAN)
yy <- as.matrix(HAKUSAN[, c(1,2,4)])
nc <- dim(yy)[1]
n <- seq(1, nc, by = 2) 
y <- yy[n, ]
z <- marfit(y, lag = 20)

marspc(z$arcoef, v = z$v)
}

\keyword{ts}
