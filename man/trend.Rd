\name{trend}
\alias{trend}

\title{Trend Estimation}
\description{
  Estimate the trend by state space model.
}
\usage{
trend(y, trend.order = 1, tau2.ini = NULL, delta, plot = TRUE, \dots)
}
\arguments{
  \item{y}{a univariate time series.}
  \item{trend.order}{trend order.}
  \item{tau2.ini}{initial estimate of variance of the system noise \eqn{\tau^2}.
    If \code{tau2.ini} = \code{NULL}, the most suitable value is chosen in
    \eqn{\tau^2 = 2^{-k}}.}
  \item{delta}{search width (for \code{tau2.ini} is specified (not \code{NULL}))
    .}
  \item{plot}{logical. If \code{TRUE} (default), '\code{trend}' and
    '\code{residual}' are plotted.}
  \item{\dots}{further arguments to be passed to \link{plot.trend}.}}
\details{
  The trend model can be represented by a state space model

    \deqn{x_n = Fx_{n-1} + Gv_n,}{x(n) = F * x(n-1) + G * v(n),}
    \deqn{y_n = Hx_n + w_n,}{y(n) = H * x(n) + w(n),}

  where \eqn{F}, \eqn{G} and \eqn{H} are matrices with appropriate dimensions.
  We assume that \eqn{v_n}{v(n)} and \eqn{w_n}{w(n)} are white noises that have
  the normal distributions \eqn{N(0,\tau^2)} and \eqn{N(0, \sigma^2)},
  respectively.
}
\value{
  An object of class \code{"trend"}, which is a list with the following
  elements:
  \item{trend}{trend component.}
  \item{residual}{residuals.}
  \item{tau2}{variance of the system noise \eqn{\tau^2}.}
  \item{sigma2}{variance of the observational noise \eqn{\sigma^2}.}
  \item{llkhood}{log-likelihood of the model.}
  \item{aic}{AIC.}
}
\references{
  Kitagawa, G. (2010)
  \emph{Introduction to Time Series Modeling}. Chapman & Hall/CRC.
}
\examples{
# The daily maximum temperatures for Tokyo
data(Temperature)
trend(Temperature, trend.order = 1, tau2.ini = 0.223, delta = 0.001)

trend(Temperature, trend.order = 2)
}

\keyword{ts}
