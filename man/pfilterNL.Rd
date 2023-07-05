\name{pfilterNL}
\alias{pfilterNL}

\title{Particle Filtering and Smoothing for Nonlinear State-Space Model}
\description{
  Trend estimation by particle filter and smoother via nonlinear state-space model.
}
\usage{
pfilterNL(y, m = 10000, lag = 20, sigma2, tau2, xrange = NULL, seed = NULL,
          plot = TRUE, \dots)
}
\arguments{
  \item{y}{univariate time series.}
  \item{m}{number of particles.}
  \item{lag}{lag length for fixed-lag smoothing.}
  \item{sigma2}{observation noise variance.}
  \item{tau2}{system noise variance.}
  \item{xrange}{specify the lower and upper bounds of the distribution's range.} 
  \item{seed}{arbitrary positive integer to generate a sequence of uniform
    random numbers. The default seed is based on the current time.}
  \item{plot}{logical. If \code{TRUE} (default), marginal smoothed distribution
    is plotted.}
  \item{\dots}{graphical arguments passed to the \code{plot} method.}
}
\details{
  This function performs particle filtering and smoothing for the following nonlinear state-space model;

  \tabular{ll}{
    \eqn{x_n = \frac{1}{2} x_{n-1} + \frac{25 x_{n-1}}{x_{n-1}^2 + 1} + 8cos(1.2n) + v_n,}{%
         x(n) = x(n-1) / 2 + (25x(n-1)) / (x(n-1)**2 + 1) + 8cos(1.2n) + v(n),} \tab (system model) \cr
    \eqn{y_n = \frac{x_n^2}{10} + w_n,}{y(n) = x(n)**2 / 10 + w(n),} \tab (observation model)
  }

where \eqn{y_n}{y(n)} is a time series, \eqn{x_n}{x(n)} is the state vector.
The system noise \eqn{v_n}{v(n)} and the observation noise \eqn{w_n}{w(n)} are assumed to be white noises
 which follow a Gaussian distribution and \eqn{v_0}{v(0)} ~ \eqn{N(0, 5)}.


The algorithm of the particle filtering and smoothing are presented in Kitagawa (2020).
For more details, please refer to Kitagawa (1996) and Doucet et al. (2001).
}
\seealso{
\code{\link{pfilter}} performs particle filtering and smoothing for linear
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
data(NLmodel)
x <- NLmodel[, 2]
pfilterNL(x, m = 100000, lag = 20 , sigma2 = 10.0, tau2 = 1.0,
          xrange = c(-20, 20), seed = 2019071117)
}

\keyword{ts}
