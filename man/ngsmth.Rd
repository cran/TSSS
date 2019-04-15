\name{ngsmth}
\alias{ngsmth}

\title{Non-Gaussian Smoothing}
\description{
  Trend estimation by non-Gaussian smoothing.
}
\usage{
ngsmth(y, noisev = 2, tau2, bv = 1.0, noisew = 1, sigma2, bw = 1.0, 
       initd = 1, k = 200, plot = TRUE, \dots)
}
\arguments{
  \item{y}{a univariate time series.}
  \item{noisev}{type of system noise density.
    \tabular{rl}{
      1 : \tab Gaussian (normal) \cr
      2 : \tab Pearson family \cr
      3 : \tab two-sided exponential}
  }
  \item{tau2}{variance of dispersion of system noise.}
  \item{bv}{shape parameter of system noise (for \code{noisev} = 2).}
  \item{noisew}{type of observation noise density
    \tabular{rl}{
      1 : \tab Gaussian (normal) \cr
      2 : \tab Pearson family \cr
      3 : \tab two-sided exponential \cr
      4 : \tab double exponential}
  }
  \item{sigma2}{variance of dispersion of observation noise.}
  \item{bw}{shape parameter of observation noise (for \code{noisew} = 2).}
  \item{initd}{type of density function.
    \tabular{rl}{
      1 : \tab Gaussian (normal) \cr
      2 : \tab uniform \cr
      3 : \tab two-sided exponential}
  }
  \item{k}{number of intervals.}
  \item{plot}{logical. If \code{TRUE} (default), '\code{trend}' and '\code{smt}'
              are plotted.}
  \item{\dots}{further arguments to be passed to \code{\link{plot.ngsmth}}.}
}
\details{
  Consider a one-dimensional state space model
    \deqn{x_n = x_{n-1} + v_n,}{x(n) = x(n-1) + v(n),}
    \deqn{y_n = x_n + w_n,}{y(n) = x(n) + w(n),}
  where the observation noise \eqn{w_n}{w(n)} is assumed to be Gaussian
  distributed and the system noise \eqn{v_n}{v(n)} is assumed to be distributed
  as the Pearson system
    \deqn{q(v_n) = c/(\tau^2 + v_n^2)^b}{q(v(n)) = c/{\tau^2 + v(n)^2}^b}
  with \eqn{\frac{1}{2} < b < \infty}{1/2 < b < \infty} and
  \eqn{c = \tau^{2b-1}\Gamma(b) \bigm/ \Gamma(\frac{1}{2})\Gamma(b-\frac{1}{2})}{%
  c = \tau^(2b-1) \Gamma(b) / {\Gamma(1/2) \Gamma(b-1/2)}}.

  This broad family of distributions includes the Cauchy distribution
  (\eqn{b = 1}) and \eqn{t}-distribution (\eqn{b = (k+1)/2}).
}
\value{
  An object of class \code{"ngsmth"}. It contains the following components:
  \item{trend}{trend.}
  \item{smt}{smoothed density.}
}
\references{
  Kitagawa, G. (2010)
  \emph{Introduction to Time Series Modeling}. Chapman & Hall/CRC.

  Kitagawa, G. and Gersch, W. (1996)
  \emph{Smoothness Priors Analysis of Time Series}. Lecture Notes in Statistics,
  No.116, Springer-Verlag.
}
\examples{
## trend model
x <- rep(0, 400)
x[101:200] <- 1
x[201:300] <- -1
y <- x + rnorm(400, mean = 0, sd = 1.0)

# system noise density : Gaussian (normal)
s1 <- ngsmth(y, noisev = 1, tau2 = 1.4e-02, noisew = 2, sigma2 = 1.048)

plot(s1, "smt", theta = 20, phi = 60, expand = 0.3)

# system noise density : Pearson family
s2 <- ngsmth(y, noisev = 2, tau2 = 2.11e-10, bv = 0.6, noisew = 2,
             sigma2 = 1.042)

plot(s2, "smt", theta = 25, phi = 30, expand = 0.25)

## seismic data
data(MYE1F)
n <- length(MYE1F)
yy <- rep(0, n)
for (i in 2:n) yy[i] <- MYE1F[i] - 0.5 * MYE1F[i-1]
m <- seq(1, n, by = 2)
y <- yy[m]
z <- tvvar(y, trend.order = 2, tau2.ini = 4.909e-02, delta = 1.0e-06)

# system noise density : Gaussian (normal)
s3 <- ngsmth(z$sm, noisev = 1, tau2 = z$tau2, noisew = 2, sigma2 = pi*pi/6,
             k = 190)

plot(s3, "smt", phi = 50, expand = 0.5, col = 8)
}



\keyword{ts}
