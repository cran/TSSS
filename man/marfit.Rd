\name{marfit}
\alias{marfit}

\title{Yule-Walker Method of Fitting Multivariate AR Model}
\description{
  Fit a multivariate AR model by the Yule-Walker method.
}
\usage{
marfit(y, lag = NULL)
}
\arguments{
  \item{y}{a multivariate time series.}
  \item{lag}{highest order of fitted AR models. Default is \eqn{2 \sqrt{n}}{%
    2*sqrt(n)}, where \eqn{n} is the length of the time series \code{y}.}
}
\value{
  An object of class \code{"maryule"}, which is a list with the following
  components:
  \item{maice.order}{order of minimum AIC.}
  \item{aic}{AIC's of the AR models with order \eqn{0,\dots,}\code{lag}.}
  \item{v}{innovation covariance matrix of the AIC best model.}
  \item{arcoef}{AR coefficients of the AIC best model.}
}
\references{
  Kitagawa, G. (2020)
  \emph{Introduction to Time Series Modeling with Applications in R}.
  Chapman & Hall/CRC.
}
\examples{
# Yaw rate, rolling, pitching and rudder angle of a ship
data(HAKUSAN)
yy <- as.matrix(HAKUSAN[, c(1,2,4)])   # Yaw rate, Pitching, Rudder angle
nc <- dim(yy)[1]
n <- seq(1, nc, by = 2) 
y <- yy[n, ]
marfit(y, 20)
}

\keyword{ts}
