\name{armachar}
\alias{armachar}

\title{Calculate Characteristics of Scalar ARMA Model}
\description{
  Calculate impulse response function, autocovariance function, autocorrelation
  function and characteristic roots of given scalar ARMA model.
}
\usage{
armachar(arcoef = NULL, macoef = NULL, v, lag = 50, nf = 200, plot = TRUE, \dots)
}
\arguments{
  \item{arcoef}{AR coefficients.}
  \item{macoef}{MA coefficients.}
  \item{v}{innovation variance.}
  \item{lag}{maximum lag of autocovariance function.}
  \item{nf}{number of frequencies in evaluating spectrum.}
  \item{plot}{logical. If \code{TRUE} (default), impulse response function,
   autocovariance, power spectrum, PARCOR and characteristic roots are plotted.}
  \item{\dots}{graphical arguments passed to the \code{plot} method.}
}
\details{
  The ARMA model is given by

  \deqn{y_t - a_1y_{t-1} - \dots - a_py_{t-p} = u_t - b_1u_{t-1} - \dots - b_qu_{t-q},}{%
        y(t) - a(1)y(t-1) - \dots - a(p)y(t-p) = u(t) - b(1)u(t-1) - \dots - b(q)u(t-q),}

  where \eqn{p} is AR order, \eqn{q} is MA order and \eqn{u_t}{u(t)} is a zero
  mean white noise.

  Characteristic roots of AR / MA operator is a list with the following components:
  \itemize{
    \item re: real part \eqn{R}
    \item im: imaginary part \eqn{I}
    \item amp: \eqn{\sqrt{R^2+I^2}}{sqrt(R^2+I^2)}
    \item atan: \eqn{\arctan(I/R)}{atan(I/R)}
    \item degree
  }
}
\value{
  An object of class \code{"arma"} which has a \code{plot} method. This is a
  list with components:
  \item{impuls}{impulse response function.}
  \item{acov}{autocovariance function.}
  \item{parcor}{PARCOR.}
  \item{spec}{power spectrum.}
  \item{croot.ar}{characteristic roots of AR operator. See Details.}
  \item{croot.ma}{characteristic roots of MA operator. See Details.}
}
\references{
  Kitagawa, G. (2020)
  \emph{Introduction to Time Series Modeling with Applications in R}.
  Chapman & Hall/CRC.
}
\examples{
# AR model : y(n) = a(1)*y(n-1) + a(2)*y(n-2) + v(n)
a <- c(0.9 * sqrt(3), -0.81)
armachar(arcoef = a, v = 1.0, lag = 20)

# MA model : y(n) = v(n) - b(1)*v(n-1) - b(2)*v(n-2)
b <- c(0.9 * sqrt(2), -0.81)
armachar(macoef = b, v = 1.0, lag = 20)

# ARMA model :  y(n) = a(1)*y(n-1) + a(2)*y(n-2)
#                      + v(n) - b(1)*v(n-1) - b(2)*v(n-2)
armachar(arcoef = a, macoef = b, v = 1.0, lag = 20)
}

\keyword{ts}
