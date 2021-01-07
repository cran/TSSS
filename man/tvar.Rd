\name{tvar}
\alias{tvar}

\title{Time Varying Coefficients AR Model}
\description{
  Estimate time varying coefficients AR model.
}
\usage{
tvar(y, trend.order = 2, ar.order = 2, span, outlier = NULL, tau2.ini = NULL,
     delta, plot = TRUE)
}
\arguments{
  \item{y}{a univariate time series.}
  \item{trend.order}{trend order (1 or 2).}
  \item{ar.order}{AR order.}
  \item{span}{local stationary span.}
  \item{outlier}{positions of outliers.}
  \item{tau2.ini}{initial estimate of variance of the system noise \eqn{\tau^2}.
    If \code{tau2.ini} = \code{NULL}, the most suitable value is chosen in
    \eqn{\tau^2 = 2^{-k}}.}
  \item{delta}{search width.} 
  \item{plot}{logical. If \code{TRUE} (default), PARCOR is plotted.}
}
\details{
  The time-varying coefficients AR model is given by

    \deqn{y_t = a_{1,t}y_{t-1} + \ldots + a_{p,t}y_{t-p} + u_t}{%
          y(t) = a(1,t)y(t-1) + ... + a(p,t)y(t-p) + u(t)}

    where \eqn{a_{i,t}}{a(i,t)} is \eqn{i}-lag AR coefficient at time \eqn{t}
    and \eqn{u_t}{u(t)} is a zero mean white noise.

  The time-varying spectrum can be plotted using AR coefficient \code{arcoef}
  and variance of the observational noise \code{sigma2} by \code{\link{tvspc}}.
}
\value{
  \item{arcoef}{time varying AR coefficients.}
  \item{sigma2}{variance of the observational noise \eqn{\sigma^2}.}
  \item{tau2}{variance of the system noise \eqn{\tau^2}.}
  \item{llkhood}{log-likelihood of the model.}
  \item{aic}{AIC.}
  \item{parcor}{PARCOR.}
}
\references{
  Kitagawa, G. (2020)
  \emph{Introduction to Time Series Modeling with Applications in R}.
  Chapman & Hall/CRC.

  Kitagawa, G. and Gersch, W. (1996)
  \emph{Smoothness Priors Analysis of Time Series}.
  Lecture Notes in Statistics, No.116, Springer-Verlag.

  Kitagawa, G. and Gersch, W. (1985)
  \emph{A smoothness priors time varying AR coefficient modeling of
  nonstationary time series}. IEEE trans. on Automatic Control, AC-30, 48-56.
}

\seealso{
\code{\link{tvspc}}, \code{\link{plot.tvspc}}
}

\examples{
# seismic data
data(MYE1F)
z <- tvar(MYE1F, trend.order = 2, ar.order = 8, span = 20,
          outlier = c(630, 1026), tau2.ini = 6.6e-06, delta = 1.0e-06)
z

spec <- tvspc(z$arcoef, z$sigma2)
plot(spec)
}

\keyword{ts}
