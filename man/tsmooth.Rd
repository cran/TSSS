\name{tsmooth}
\alias{tsmooth}

\title{Prediction and Interpolation of Time Series}
\description{
Predict and interpolate time series based on state space model by Kalman filter.
}
\usage{
tsmooth(y, f, g, h, q, r, x0 = NULL, v0 = NULL, filter.end = NULL,
        predict.end = NULL, minmax = c(-1.0e+30, 1.0e+30), missed = NULL,
        np = NULL, plot = TRUE, \dots)
}
\arguments{
  \item{y}{a univariate time series \eqn{y_n}{y(n)}.}
  \item{f}{state transition matrix \eqn{F_n}{F(n)}.}
  \item{g}{matrix \eqn{G_n}{G(n)}.}
  \item{h}{matrix \eqn{H_n}{H(n)}.}
  \item{q}{system noise variance \eqn{Q_n}{Q(n)}.}
  \item{r}{observational noise variance \eqn{R}.}
  \item{x0}{initial state vector \eqn{X(0\mid0)}{X(0|0)}.}
  \item{v0}{initial state covariance matrix \eqn{V(0\mid0)}{V(0|0)}.}
  \item{filter.end}{end point of filtering.}
  \item{predict.end}{end point of prediction.}
  \item{minmax}{lower and upper limits of observations.}
  \item{missed}{start position of missed intervals.}
  \item{np}{number of missed observations.}
  \item{plot}{logical. If \code{TRUE} (default), '\code{mean.smooth}' and
    '\code{esterr}' are plotted.}
  \item{\dots}{further arguments to be passed to \code{\link{plot.smooth}}.}}
\details{
  The linear Gaussian state space model is

    \deqn{x_n = F_n x_{n-1} + G_n v_n,}{x(n) = F(n)x(n-1) + G(n)v(n),}
    \deqn{y_n = H_n x_n + w_n,}{y(n) = H(n)x(n) + w(n),}

  where \eqn{y_n}{y(n)} is a univariate time series, \eqn{x_n}{x(n)} is an
    \eqn{m}-dimensional state vector.
  
  \eqn{F_n}{F(n)}, \eqn{G_n}{G(n)} and \eqn{H_n}{H(n)} are \eqn{m \times m}{m * m},
  \eqn{m \times k}{m * k} matrices and a vector of length \eqn{m} , respectively.
  \eqn{Q_n}{Q(n)} is \eqn{k \times k}{k * k} matrix and \eqn{R_n}{R(n)} is a scalar.
  \eqn{v_n}{v(n)} is system noise and \eqn{w_n}{w(n)} is observation noise,
  where we assume that \eqn{E(v_n, w_n) = 0}{E(v(n), w(n)) = 0},
  \eqn{v_n \sim N(0, Q_n)}{v(n) ~ N(0, Q(n))} and
  \eqn{w_n \sim N(0, R_n)}{w(n) ~ N(0, R(n))}. User should give all the matrices
  of a state space model and its parameters. In current version, \eqn{F_n}{F(n)},
  \eqn{G_n}{G(n)}, \eqn{H_n}{H(n)}, \eqn{Q_n}{Q(n)}, \eqn{R_n}{R(n)} should be
  time invariant.
}
\value{
  An object of class \code{"smooth"}. It contains the following components:
  \item{mean.smooth}{mean vectors of the smoother.}
  \item{cov.smooth}{variance of the smoother.}
  \item{esterr}{estimation error.}
  \item{llkhood}{log-likelihood.}
  \item{aic}{AIC.}
}
\references{
  Kitagawa, G. (2010)
  \emph{Introduction to Time Series Modeling}. Chapman & Hall/CRC.

  Kitagawa, G. and Gersch, W. (1996)
  \emph{Smoothness Priors Analysis of Time Series}.
  Lecture Notes in Statistics, No.116, Springer-Verlag.
}
\examples{
## Example of prediction (AR model : m=15, k=1)
data(BLSALLFOOD)
BLS120 <- BLSALLFOOD[1:120]
z1 <- arfit(BLS120, plot = FALSE)
tau2 <- z1$sigma2
arcoef <- z1$arcoef

# in case m = 15
m1 <- z1$maice.order
f <- matrix(0.0e0, m1, m1)
f[1, ] <- arcoef[1:m1]
if (m1 != 1)
  for (i in 2:m1) f[i, i-1] <- 1
g <- c(1, rep(0.0e0, m1-1))
h <- c(1, rep(0.0e0, m1-1))
q <- tau2[m1+1]
r <- 0.0e0
x0 <- rep(0.0e0, m1)
v0 <- NULL

s1 <- tsmooth(BLS120, f, g, h, q, r, x0, v0, filter.end = 120, predict.end = 156)
s1

plot(s1, BLSALLFOOD)

## Example of interpolation of missing values  (AR model : m=15, k=1)
z2 <- arfit(BLSALLFOOD, plot = FALSE)
tau2 <- z2$sigma2
arcoef <- z2$arcoef

# in case m2 = 15
m2 <- z2$maice.order
f <- matrix(0.0e0, m2, m2)
f[1, ] <- arcoef[1:m2]
if (m2 != 1)
  for (i in 2:m2) f[i, i-1] <- 1
g <- c(1, rep(0.0e0, m2-1))
h <- c(1, rep(0.0e0, m2-1))
q <- tau2[m2+1]
r <- 0.0e0
x0 <- rep(0.0e0, m2)
v0 <- NULL

tsmooth(BLSALLFOOD, f, g, h, q, r, x0, v0, missed = c(41, 101), np = c(30, 20))
}

\keyword{ts}
