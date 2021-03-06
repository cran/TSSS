\name{lsar}
\alias{lsar}

\title{Decomposition of Time Interval to Stationary Subintervals}
\description{
  Decompose time series to stationary subintervals and estimate local spectrum.
}
\usage{
lsar(y, max.arorder = 20, ns0, plot = TRUE, \dots)
}
\arguments{
  \item{y}{a univariate time series.}
  \item{max.arorder}{highest order of AR model.}
  \item{ns0}{length of basic local span.}
  \item{plot}{logical. If \code{TRUE} (default), local spectra are plotted.}
  \item{\dots}{graphical arguments passed to the \code{plot} method.}
}
\value{
  An object of class \code{"lsar"} which has a \code{plot} method. This is a
  list with the following components:
  \item{model}{1: pooled model is accepted.\cr 2: switched model is accepted.} 
  \item{ns}{number of observations of local span.}
  \item{span}{start points and end points of local spans.}
  \item{nf}{number of frequencies in computing local power spectrum.}
  \item{ms}{order of switched model.}
  \item{sds}{innovation variance of switched model.}
  \item{aics}{AIC of switched model.}
  \item{mp}{order of pooled model.}
  \item{sdp}{innovation variance of pooled model.}
  \item{aics}{AIC of pooled model.}
  \item{spec}{local spectrum.}
  \item{tsname}{the name of the univariate time series \code{y}.}
}
\references{
  Kitagawa, G. (2020)
  \emph{Introduction to Time Series Modeling with Applications in R}.
  Chapman & Hall/CRC.
}
\examples{
# seismic data
data(MYE1F)
lsar(MYE1F, max.arorder = 10, ns0 = 100)
}

\keyword{ts}
