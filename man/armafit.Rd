\name{armafit}
\alias{armafit}

\title{Scalar ARMA Model Fitting}
\description{
Fit a scalar ARMA model by maximum likelihood method.
}
\usage{
armafit(y, ar.order, ar = NULL, ma.order, ma = NULL)
}
\arguments{
  \item{y}{a univariate time series.}
  \item{ar.order}{AR order.}
  \item{ar}{initial AR coefficients. If \code{NULL} (default), use default
    initial values.}
  \item{ma.order}{MA order.}
  \item{ma}{initial MA coefficients. If \code{NULL} (default), use default
    initial values.}
}
\value{
  \item{sigma2}{innovation variance.}
  \item{llkhood}{log-likelihood of the model.}
  \item{aic}{AIC of the model.}
  \item{arcoef}{AR coefficients.}
  \item{macoef}{MA coefficients.}
}
\references{
  Kitagawa, G. (2020)
  \emph{Introduction to Time Series Modeling with Applications in R}.
  Chapman & Hall/CRC.
}
\examples{
# Sunspot number data
data(Sunspot)
y <- log10(Sunspot)
z <- armafit(y, ar.order = 3, ma.order = 3)
z

armachar(arcoef = z$arcoef, macoef = z$macoef, v = z$sigma2, lag = 20)
}

\keyword{ts}
