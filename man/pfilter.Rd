\name{pfilter}
\alias{pfilter}
\alias{plot.pfilter}

\title{Particle Filtering and Smoothing}
\description{
  Trend estimation by particle filter and smoother.
}
\usage{
pfilter(y, m = 10000, model = 0, lag = 20, initd = 0, sigma2, tau2,
        alpha = 0.99, bigtau2 = NULL, init.sigma2 = 1, xrange = NULL,
        seed = NULL, plot = TRUE, \dots)
}
\arguments{
  \item{y}{univariate time series.}
  \item{m}{number of particles.}
  \item{model}{model for the system noise.
    \tabular{rl}{
      0: \tab normal distribution \cr
      1: \tab Cauchy distribution \cr
      2: \tab Gaussian mixture distribution \cr
           \tab \eqn{\alpha} \eqn{N}(0, \eqn{\tau^2}) + (1 - \eqn{\alpha}) \eqn{N}(0, \eqn{T^2}), \cr
           \tab where \eqn{N} is the normal density.}}
  \item{lag}{lag length for fixed-lag smoothing.}
  \item{initd}{type of initial state distribution.
    \tabular{rl}{
      0: \tab normal distribution \cr
      1: \tab uniform distribution \cr
      2: \tab Cauchy distribution \cr
      3: \tab fixed point (default value = 0)}}
  \item{sigma2}{observation noise variance \eqn{\sigma^2}.}
  \item{tau2}{system noise variance \eqn{\tau^2} for \code{model} = 0 or dispersion parameter
    for \code{model} = 1.}
  \item{alpha}{mixture weight \eqn{\alpha}. (valid for \code{model} = 2)}
  \item{bigtau2}{variance of the second component \eqn{T^2}. (valid for \code{model} = 2)}
  \item{init.sigma2}{variance for \code{initd} = 0 or dispersion parameter of
    initial state distribution for \code{initd} = 2.}
  \item{xrange}{specify the lower and upper bounds of the distribution's range.} 
  \item{seed}{arbitrary positive integer to generate a sequence of uniform
    random numbers. The default seed is based on the current time.}
  \item{plot}{logical. If \code{TRUE} (default), marginal smoothed distribution
              is plotted.}
  \item{\dots}{graphical arguments passed to the \code{plot} method.}
}
\details{
This function performs particle filtering and smoothing for the first order trend model;

  \tabular{ll}{
    \eqn{x_n = x_{n-1} + v_n,}{x(n) = x(n-1) + v(n),} \tab (system model) \cr 
    \eqn{y_n = x_n + w_n,}{y(n) = x(n) + w(n),} \tab (observation model)} 

where \eqn{y_n}{y(n)} is a time series, \eqn{x_n}{x(n)} is the state vector.
The system noise \eqn{v_n}{v(n)} and the observation noise \eqn{w_n}{w(n)} are assumed to be white noises
 which follow a Gaussian distribution or a Cauchy distribution, and non-Gaussian distribution, respectively.


The algorithm of the particle filter and smoother are presented in Kitagawa (2020).
For more details, please refer to Kitagawa (1996) and Doucet et al. (2001).
}
\seealso{
\code{\link{pfilterNL}} performs particle filtering and smoothing for nonlinear
 non-Gaussian state-space model.
}
\value{
  An object of class \code{"pfilter"} which has a \code{plot} method. This is a
  list with the following components:
  \item{llkhood}{log-likelihood.}
  \item{smooth.dist}{marginal smoothed distribution of the trend \eqn{T(i,j)} 
      \eqn{(i = 1,...,n, j = 1,...,7)}, where \eqn{n} is the length of \code{y}.
    \tabular{rll}{
      \tab j = 4: \tab 50\% point \cr
      \tab j = 3, 5: \tab 1-sigma points (15.87\% and 84.14\% points) \cr
      \tab j = 2, 6: \tab 2-sigma points (2.27\% and 97.73\% points) \cr
      \tab j = 1, 7: \tab 3-sigma points (0.13\% and 99.87\% points)}}
}
\references{
  Kitagawa, G. (1996)
  \emph{Monte Carlo filter and smoother for non-Gaussian nonlinear state space models},
 J. of Comp. and Graph. Statist., 5, 1-25.

  Doucet, A., de Freitas, N. and Gordon, N. (2001)
  \emph{Sequential Monte Carlo Methods in Practice}, Springer, New York.

  Kitagawa, G. (2020)
  \emph{Introduction to Time Series Modeling with Applications in R}. Chapman & Hall/CRC.
}
\examples{
data(PfilterSample)
y <- PfilterSample

\dontrun{
pfilter(y, m = 100000, model = 0, lag = 20, initd = 0, sigma2 = 1.048,
        tau2 = 1.4e-2, xrange = c(-4, 4), seed = 2019071117)

pfilter(y, m = 100000, model = 1, lag = 20 , initd = 0, sigma2 = 1.045,
        tau2 = 3.53e-5, xrange = c(-4, 4), seed = 2019071117)

pfilter(y, m = 100000, model = 2, lag = 20 , initd = 0, sigma2 = 1.03,
        tau2 = 0.00013, alpha = 0.991, xrange = c(-4, 4), seed = 2019071117)
}
}

\keyword{ts}
