\name{ngsim}
\alias{ngsim}

\title{Simulation by Non-Gaussian State Space Model}
\description{
  Simulation by non-Gaussian state space model.
}
\usage{
ngsim(n = 200, trend = NULL, seasonal.order = 0, seasonal = NULL, arcoef = NULL, 
      ar = NULL, noisew = 1, wminmax = c(-1, 1), paramw = NULL, noisev = 1, 
      vminmax = c(-1, 1), paramv = NULL, seed = NULL, plot = TRUE, \dots)
}
\arguments{
  \item{n}{the number of simulated data.}
  \item{trend}{initial values of trend component of length at most 2.}
  \item{seasonal.order}{seasonal order. (0 or 1)} 
  \item{seasonal}{if \code{seasonal.order} > 0, initial values of seasonal
    component of length \eqn{p-1}, where \eqn{p} is the number of season in one
    period.}
  \item{arcoef}{AR coefficients.}
  \item{ar}{initial values of AR component.}
  \item{noisew}{type of the observational noise.
    \tabular{rl}{
      -1 : \tab Cauchy random number (without an inverse function) \cr
      -2 : \tab exponential distribution (without an inverse function)\cr
      -3 : \tab double exponential distribution (without an inverse function)\cr
       0 : \tab double exponential distribution (+ Euler's constant)\cr
       1 : \tab normal distribution,\cr
       2 : \tab Pearson distribution,\cr
       3 : \tab double exponential distribution}
  }
  \item{wminmax}{lower and upper bound of observational noise.}
  \item{paramw}{parameter of the observational noise density.
    \tabular{rl}{
      noisew = 1 : \tab variance\cr
      noisew = 2 : \tab dispersion parameter (tau square), shape parameter}
  }
  \item{noisev}{type of the system noise.
    \tabular{rl}{
      -1 : \tab Cauchy random number (without an inverse function)\cr
      -2 : \tab exponential distribution (without an inverse function)\cr
      -3 : \tab double exponential distribution (without an inverse function)\cr
       0 : \tab double exponential distribution (+ Euler's constant)\cr
       1 : \tab normal distribution\cr
       2 : \tab Pearson distribution\cr
       3 : \tab double exponential distribution}
  }
  \item{vminmax}{lower and upper bound of system noise.}
  \item{paramv}{parameter of the system noise density.
    \tabular{rl}{
    noisev = 1 : \tab variance\cr
    noisev = 2 : \tab dispersion parameter (tau square), shape parameter}
  }
  \item{seed}{arbitrary positive integer to generate a sequence of uniform
    random numbers. The default seed is based on the current time.}
  \item{plot}{logical. If \code{TRUE} (default), simulated data are plotted.}
  \item{\dots}{further arguments to be passed to \code{\link{plot.simulate}}.}
}
\value{
  An object of class \code{"simulate"}, giving simulated data of non-Gaussian
  state space model.
}
\references{
  Kitagawa, G. (2010)
  \emph{Introduction to Time Series Modeling}. Chapman & Hall/CRC.
}
\examples{
ar1 <- ngsim(n = 400, arcoef = 0.95, noisew = 1, paramw = 1, noisev = 1,
             paramv = 1, seed = 555)
plot(ar1, use = c(201, 400))

ar2 <- ngsim(n = 400, arcoef = c(1.3, -0.8), noisew = 1, paramw = 1, noisev = 1,
             paramv = 1, seed = 555)
plot(ar2, use = c(201, 400))


}

\keyword{ts}
