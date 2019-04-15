\name{marspc}
\alias{marspc}

\title{Cross Spectra and Power Contribution}
\description{
  Compute cross spectra and power contribution.
}
\usage{
marspc(arcoef, v, plot = TRUE, \dots)
}
\arguments{
  \item{arcoef}{AR coefficient matrices.}
  \item{v}{innovation variance matrix.}
  \item{plot}{logical. If \code{TRUE} (default), cross spectra and power
    contribution are plotted.}
  \item{\dots}{further arguments to be passed to \code{plot.marspc}.}
}
\value{
  An object of class \code{"marspc"}, which is a list with the following
  elements:
  \item{spec}{cross spectra.}
  \item{amp}{amplitude spectra.}
  \item{phase}{Phase spectra.}
  \item{coh}{simple coherency.}
  \item{power}{power contribution.}
  \item{rpower}{relative power contribution.}
}
\references{
  Kitagawa, G. (2010)
  \emph{Introduction to Time Series Modeling}. Chapman & Hall/CRC.
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
