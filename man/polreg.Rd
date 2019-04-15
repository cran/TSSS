\name{polreg}
\alias{polreg}

\title{Polynomial Regression Model}
\description{
  Estimate the trend using the AIC best polynomial regression model.
}
\usage{
polreg(y, order, plot = TRUE, \dots)
}
\arguments{
  \item{y}{a univariate time series.}
  \item{order}{order of polynomial regression.}
  \item{plot}{logical. If \code{TRUE} (default), '\code{y}' and '\code{trend}'
    are plotted.}
  \item{\dots}{further arguments to be passed to \link{plot.polreg}.}
}
\value{
  An object of class \code{"polreg"}, which is a list with the following
  elements:
  \item{order.maice}{MAICE (minimum AIC estimate) order.}
  \item{sigma2}{residual variance of the model with order \eqn{M}.
    (\eqn{0 \leq M \leq}{0 <= M <=} \code{order} \eqn{+ 1})}
  \item{aic}{AIC of the model with order \eqn{M}.
    (\eqn{0 \leq M \leq}{0 <= M <=} \code{order} \eqn{+ 1})}
  \item{daic}{AIC - minimum AIC.}
  \item{coef}{regression coefficients \eqn{A(I,M)} with order \eqn{M}.

    (\eqn{1 \leq M \leq}{1 <= M <=} \code{order} \eqn{+ 1},
    \eqn{1 \leq I \leq M}{1 <= I <= M})}
  \item{trend}{trend component.}
}
\references{
  Kitagawa, G. (2010)
  \emph{Introduction to Time Series Modeling}. Chapman & Hall/CRC.
}
\examples{
# The daily maximum temperatures for Tokyo
data(Temperature)
polreg(Temperature, order = 7)

# Wholesale hardware data
data(WHARD)
y <- log10(WHARD)
polreg(y, order = 15)
}

\keyword{ts}
