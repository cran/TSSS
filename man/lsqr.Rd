\name{lsqr}
\alias{lsqr}

\title{The Least Squares Method via Householder Transformation}
\description{
  Compute Regression coefficients of the model with minimum AIC.
}
\usage{
lsqr(y, lag = 10, plot = TRUE, \dots)
}
\arguments{
  \item{y}{a univariate time series.}
  \item{lag}{number of sine and cosine terms.}
  \item{plot}{logical. If \code{TRUE} (default), original data and fitted
    trigonometric polynomial are plotted.}
  \item{\dots}{further arguments to be passed to \code{\link{plot.lsqr}}.}
}
\value{
  An object of class \code{"lsqr"}, which is a list with the following
  elements:
  \item{aic}{
    AIC's of the model with order \eqn{0,\dots,k ( = 2}\code{lag}\eqn{ + 1)}.}
  \item{sigma2}{residual variance of the model with order \eqn{0,\dots,k}.}
  \item{maice.order}{order of minimum AIC.}
  \item{regress}{regression coefficients of the model.}
  \item{tripoly}{trigonometric polynomial.}
}
\references{
  Kitagawa, G. (2010)
  \emph{Introduction to Time Series Modeling}. Chapman & Hall/CRC.
}
\examples{
# The daily maximum temperatures for Tokyo
data(Temperature)
lsqr(Temperature)
}

\keyword{ts}
