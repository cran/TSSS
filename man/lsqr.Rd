\name{lsqr}
\alias{lsqr}

\title{The Least Squares Method via Householder Transformation}
\description{
  Compute regression coefficients of the model with minimum AIC by the least
  squares method via Householder transformation.
}
\usage{
lsqr(y, lag = NULL, period = 365, plot = TRUE, \dots)
}
\arguments{
  \item{y}{a univariate time series.}
  \item{lag}{number of sine and cosine components. Default is \eqn{\sqrt{n}}{sqrt(n)},
    where \eqn{n} is the length of the time series \code{y}.}
  \item{period}{period of one cycle.}
  \item{plot}{logical. If \code{TRUE} (default), original data and fitted
    trigonometric polynomial are plotted.}
  \item{\dots}{graphical arguments passed to \code{\link{plot.lsqr}}.}
}
\value{
  An object of class \code{"lsqr"}, which is a list with the following
   components:
  \item{aic}{
    AIC's of the model with order \eqn{0,\dots,k ( = 2}\code{lag}\eqn{ + 1)}.}
  \item{sigma2}{residual variance of the model with order \eqn{0,\dots,k}.}
  \item{maice.order}{order of minimum AIC.}
  \item{regress}{regression coefficients of the model.}
  \item{tripoly}{trigonometric polynomial.}
}
\references{
  Kitagawa, G. (2020)
  \emph{Introduction to Time Series Modeling with Applications in R}.
  Chapman & Hall/CRC.
}
\examples{
# The daily maximum temperatures in Tokyo
data(Temperature)
lsqr(Temperature, lag = 10)
}

\keyword{ts}
