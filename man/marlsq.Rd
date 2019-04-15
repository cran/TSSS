\name{marlsq}
\alias{marlsq}

\title{Least Squares Method for Multivariate AR Model}
\description{
  Fit a multivariate AR model by least squares method.
}
\usage{
marlsq(y, lag = NULL)
}
\arguments{
  \item{y}{a multivariate time series.}
  \item{lag}{highest AR order. Default is \eqn{2 \sqrt{n}}{2*sqrt(n)}, where
    \eqn{n} is the length of the time series \code{y}.}
}
\value{
  An object of class \code{"marlsq"}, which is a list with the following
  elements:
  \item{maice.order}{order of the MAICE model.}
  \item{aic}{total AIC of the model.}
  \item{v}{innovation covariance matrix.}
  \item{arcoef}{AR coefficient matrices.}
}
\references{
  Kitagawa, G. (2010)
  \emph{Introduction to Time Series Modeling}. Chapman & Hall/CRC.
}
\examples{
# Yaw rate, rolling, pitching and rudder angle of a ship
data(HAKUSAN)
y <- as.matrix(HAKUSAN[, c(1,2,4)])   # Yaw rate, Rolling, Rudder angle
z <- marlsq(y)
z

marspc(z$arcoef, v = z$v)
}

\keyword{ts}
