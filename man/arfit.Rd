\name{arfit}
\alias{arfit}

\title{Univariate AR Model Fitting}
\description{
  Fit a univariate AR model by Yule-Walker method, Least squares (Householder)
  method or PARCOR method.
}
\usage{
arfit(y, lag = NULL, method = 1, plot = TRUE, \dots)
}
\arguments{
  \item{y}{a univariate time series.}
  \item{lag}{highest order of AR model. Default is \eqn{2 \sqrt{n}}{2*sqrt(n)},
    where \eqn{n} is the length of the time series \code{y}.}
  \item{method}{estimation procedure.
    \tabular{rl}{
      1 : \tab Yule-Walker method\cr
      2 : \tab Least squares (Householder) method\cr
      3 : \tab PARCOR method (Partial autoregression)\cr
      4 : \tab PARCOR method (PARCOR)\cr
      5 : \tab PARCOR method (Burg's algorithm)}
  }
  \item{plot}{logical. If \code{TRUE} (default), PARCOR, AIC and power spectrum
    are plotted.}
  \item{\dots}{further arguments to be passed to \code{plot.arfit}.}}
\value{
  An object of class \code{"arfit"}, which is a list with the following
  elements:
  \item{sigma2}{innovation variance.}
  \item{maice.order}{order of minimum AIC.}
  \item{aic}{AIC.}
  \item{arcoef}{AR coefficients of the best model.}
  \item{parcor}{PARCOR.}
  \item{spec}{power spectrum (in log scale).}
  \item{tsname}{the name of the univariate time series \code{y}.}
}
\references{
  Kitagawa, G. (2010)
  \emph{Introduction to Time Series Modeling}. Chapman & Hall/CRC.
}
\examples{
# Sun spot number data
data(Sunspot)
arfit(log10(Sunspot), 20)

# BLSALLFOOD data
data(BLSALLFOOD)
arfit(BLSALLFOOD)
}

\keyword{ts}
