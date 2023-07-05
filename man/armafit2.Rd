\name{armafit2}
\alias{armafit2}

\title{Scalar ARMA Model Fitting}
\description{
Estimate all ARMA models within the user-specified maximum order by
 maximum likelihood method.
}
\usage{
armafit2(y, ar.order, ma.order)
}
\arguments{
  \item{y}{a univariate time series.}
  \item{ar.order}{maximum AR order.}
  \item{ma.order}{maximum MA order.}
}
\value{
  \item{aicmin}{minimum AIC.}
  \item{maice.order}{AR and MA orders of minimum AIC model.}
  \item{sigma2}{innovation variance of all models.}
  \item{llkhood}{log-likelihood of all models.}
  \item{aic}{AIC of all models.}
  \item{coef}{AR and MA coefficients of all models.}
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
armafit2(y, ar.order = 5, ma.order = 5)
}

\keyword{ts}
