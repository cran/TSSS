\name{tvvar}
\alias{tvvar}
\alias{print.tvvar}
\alias{plot.tvvar}

\title{Time Varying Variance}
\description{
  Estimate time-varying variance.
}
\usage{
tvvar(y, trend.order, tau2.ini = NULL, delta, plot = TRUE, \dots)
}
\arguments{
  \item{y}{a univariate time series.}
  \item{trend.order}{trend order.}
  \item{tau2.ini}{initial estimate of variance of the system noise \eqn{\tau^2}.
    If \code{tau2.ini} = \code{NULL}, the most suitable value is chosen in
    \eqn{\tau^2 = 2^{-k}}.}
  \item{delta}{search width.}
  \item{plot}{logical. If \code{TRUE} (default), transformed data, trend and
    residuals are plotted.}
  \item{\dots}{graphical arguments passed to the \code{plot} method.}
}
\details{
  Assuming that \eqn{\sigma_{2m-1}^2 = \sigma_{2m}^2}{%
  \sigma_(2m-1)^2 = \sigma_(2m)^2}, we define a transformed time series
  \eqn{s_1,\dots,s_{N/2}}{s_1,\dots,s_(N/2)} by 

  \deqn{s_m = y_{2m-1}^2 + y_{2m}^2,}{s_m = y_(2m-1)^2 + y_(2m)^2,}

  where \eqn{y_n} is a Gaussian white noise with mean \eqn{0} and variance
  \eqn{\sigma_n^2}. \eqn{s_m} is distributed as a \eqn{\chi^2} distribution with
  \eqn{2} degrees of freedom, so the probability density function of \eqn{s_m}
  is given by

  \deqn{f(s) = \frac{1}{2\sigma^2} e^{-s/2\sigma^2}.}{%
        f(s) = e^{-s/(2\sigma^2)} / (2\sigma^2).}
 
  By further transformation

  \deqn{z_m = \log \left( \frac{s_m}{2} \right),}{z_m = log{(s_m)/2},}

  the probability density function of \eqn{z_m} is given by

  \deqn{g(z) = \frac{1}{\sigma^2} \exp{ \left\{ z-\frac{e^z}{\sigma^2} \right\} } = \exp{ \left\{ (z-\log\sigma^2) - e^{(z-\log\sigma^2)} \right\} }.}{%
        g(z) = exp{z - e^z/\sigma^2} / \sigma^2 = exp{(z-log(\sigma^2)) - e^(z-log(\sigma^2))}.}

  Therefore, the transformed time series is given by

  \deqn{z_m = \log \sigma^2 + w_m,}{z_m = log(\sigma^2) + w_m,}

  where \eqn{w_m} is a double exponential distribution with probability density
  function 

  \deqn{h(w) = \exp{\{w-e^w\}}.}{h(w) = exp{w-e^w}.}

  In the space state model 

  \deqn{z_m = t_m + w_m}

  by identifying trend components of \eqn{z_m}, the log variance of original
  time series \eqn{y_n} is obtained.
}
\value{
  An object of class \code{"tvvar"} which has a \code{plot} method. This is a
  list with the following components:
  \item{tvv}{time varying variance.}
  \item{nordata}{normalized data.}
  \item{sm}{transformed data.}
  \item{trend}{trend.}
  \item{noise}{residuals.}
  \item{tau2}{variance of the system noise.}
  \item{sigma2}{variance of the observational noise.}
  \item{llkhood}{log-likelihood of the model.}
  \item{aic}{AIC.}
  \item{tsname}{the name of the univariate time series \code{y}.}
}
\references{
  Kitagawa, G. (2020)
  \emph{Introduction to Time Series Modeling with Applications in R}.
  Chapman & Hall/CRC.

  Kitagawa, G. and Gersch, W. (1996)
  \emph{Smoothness Priors Analysis of Time Series}. Lecture Notes in Statistics,
  No.116, Springer-Verlag.

  Kitagawa, G. and Gersch, W. (1985)
  \emph{A smoothness priors time varying AR coefficient modeling of
  nonstationary time series}. IEEE trans. on Automatic Control, AC-30, 48-56.
}
\examples{
# seismic data
data(MYE1F)
tvvar(MYE1F, trend.order = 2, tau2.ini = 6.6e-06, delta = 1.0e-06)
}

\keyword{ts}
