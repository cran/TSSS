\name{season}
\alias{season}

\title{Seasonal Adjustment}
\description{
  Seasonal adjustment by state space modeling. 
}
\usage{
season(y, trend.order = 1, seasonal.order = 1, ar.order = 0, trade = FALSE,
       period = 12, tau2.ini = NULL, filter = c(1, length(y)),
       predict = length(y), arcoef.ini = NULL, log = FALSE,
       minmax = c(-1.0e+30, 1.0e+30), plot = TRUE, ...)
}
\arguments{
  \item{y}{a univariate time series.}
  \item{trend.order}{trend order.}
  \item{seasonal.order}{seasonal order.}
  \item{ar.order}{AR order.}
  \item{trade}{logical; if \code{TRUE}, the model including trading day effect
               component is considered, where \code{tsp(y)} is not \code{NULL} 
               and \code{frequency(y)} is 4 or 12.}
  \item{period}{number of seasons in one period. If the tsp attribute of 
                \code{y} is not \code{NULL}, \code{frequency(y)}.
    \tabular{rl}{
      = 12 : for monthly data\cr
      = 4 : for quarterly data
    }}
  \item{tau2.ini}{initial estimate of variance of the system noise \eqn{\tau^2},
    not equal to 1.}
  \item{filter}{a numerical vector of the form \code{c(x1,x2)} which gives start
    and end position of filtering.}
  \item{predict}{the end position of prediction (\eqn{\geq}{>=} \code{x2}).}
  \item{arcoef.ini}{initial estimate of AR coefficients (for \code{ar.order} >
    0).}
  \item{log}{logical. If \code{TRUE}, the data \code{y} is log-transformed.}
  \item{minmax}{lower and upper limits of observations.}
  \item{plot}{logical. If \code{TRUE} (default), '\code{trend}',
    '\code{seasonal}' and '\code{ar}' are plotted.}
  \item{\dots}{further arguments to be passed to \code{\link{plot.season}}.}
}
\value{
  An object of class \code{"season"}, which is a list with the following
  elements:
  \item{tau2}{variance of the system noise.}
  \item{sigma2}{variance of the observational noise.}
  \item{llkhood}{log-likelihood of the model.}
  \item{aic}{AIC of the model.}
  \item{trend}{trend component (for \code{trend.order} > 0).}
  \item{seasonal}{seasonal component (for \code{seasonal.order} > 0).}
  \item{arcoef}{AR coefficients (for \code{ar.order} > 0).}
  \item{ar}{AR component (for \code{ar.order} > 0).}
  \item{day.effect}{trading day effect (for \code{trade} = 6).}
  \item{noise}{noise component.}
  \item{cov}{covariance matrix of smoother.}
}
\references{
  Kitagawa, G. (2010)
  \emph{Introduction to Time Series Modeling}. Chapman & Hall/CRC.
}
\examples{
# BLSALLFOOD data
data(BLSALLFOOD)
season(BLSALLFOOD, trend.order = 2, seasonal.order = 1, ar.order = 2)

season(BLSALLFOOD, trend.order = 2, seasonal.order = 1, ar.order = 2,
       filter = c(1, 132))

# Wholesale hardware data
data(WHARD)
season(WHARD, trend.order = 2, seasonal.order = 1, ar.order = 0, trade = TRUE,
       log = TRUE)

season(WHARD, trend.order = 2, seasonal.order = 1, ar.order = 0, trade = TRUE,
       filter = c(1, 132), log = TRUE)
}

\keyword{ts}
