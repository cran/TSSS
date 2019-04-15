\name{boxcox}
\alias{boxcox}

\title{Box-Cox Transformation}
\description{
  Computes Box-Cox transformation and find an optimal lambda with minimum AIC.}
\usage{
boxcox(y, plot = TRUE, \dots)
}
\arguments{
  \item{y}{a univariate time series.}
  \item{plot}{logical. If \code{TRUE} (default), original data and transformed
  data with minimum AIC are plotted.}
  \item{\dots}{further arguments to be passed to \code{plot.boxcox}.}
}
\value{
  An object of class \code{"boxcox"}, which is a list with the following
  elements:
  \item{mean}{mean of original data.}
  \item{var}{variance of original data.}
  \item{aic}{AIC of the model with respect to the original data.}
  \item{llkhood}{log-likelihood of the model with respect to the original data.}
  \item{z}{transformed data.}
  \item{aic.z}{AIC of the model with respect to the transformed data.}
  \item{llkhood.z}{log-likelihood of the model with respect to the transformed
                   data.}
}
\references{
  Kitagawa, G. (2010)
  \emph{Introduction to Time Series Modeling}. Chapman & Hall/CRC.
}
\examples{
# Sun spot number data
data(Sunspot)
boxcox(Sunspot)

# Wholesale hardware data
data(WHARD)
boxcox(WHARD)
}

\keyword{ts}
