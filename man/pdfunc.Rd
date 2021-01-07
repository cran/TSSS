\name{pdfunc}
\alias{pdfunc}

\title{Probability Density Function}
\description{
  Evaluate probability density function for normal distribution, Cauchy
 distribution, Pearson distribution, exponential distribution, Chi-square
 distributions, double exponential distribution and uniform distribution.}
\usage{
pdfunc(model = "norm", mean = 0, sigma2 = 1, mu = 0, tau2 = 1, shape,
       lambda = 1, side = 1, df, xmin = 0, xmax = 1, plot = TRUE, \dots)
}
\arguments{
  \item{model}{a character string indicating the model type of probability
    density function: either \code{"norm"}, \code{"Cauchy"}, \code{"Pearson"},
    \code{"exp"}, \code{"Chi2"}, \code{"dexp"} or \code{"unif"}.}
  \item{mean}{mean. (valid for \code{"norm"})}
  \item{sigma2}{variance. (valid for \code{"norm"})}
  \item{mu}{location parameter \eqn{\mu}. (valid for \code{"Cauchy"} and
    \code{"Pearson"})}
  \item{tau2}{dispersion parameter \eqn{\tau^2}. (valid for \code{"Cauchy"} and
    \code{"Pearson"})}
  \item{shape}{shape parameter (> 0.5). (valid for \code{"Pearson"})}
  \item{lambda}{lambda \eqn{\lambda}. (valid for \code{"exp"})}
  \item{side}{1: exponential, 2: two-sided exponential.
    (valid for \code{"exp"})}
  \item{df}{degree of freedoms \eqn{k}.  (valid for \code{"Chi2"})}
  \item{xmin}{lower bound of the interval.}
  \item{xmax}{upper bound of the interval.}
  \item{plot}{logical. If \code{TRUE} (default), probability density function is
              plotted.}
  \item{\dots}{graphical arguments passed to the \code{plot} method.}
}
\value{
  An object of class \code{"pdfunc"} which has a \code{plot} method. This is a
  list with the following components:
  \item{density}{values of density function.}
  \item{interval}{lower and upper bound of interval.}
  \item{param}{parameters of model.}
}
\references{
  Kitagawa, G. (2020)
  \emph{Introduction to Time Series Modeling with Applications in R}.
  Chapman & Hall/CRC.
}
\examples{
# normal distribution
pdfunc(model = "norm", xmin = -4, xmax = 4) 

# Cauchy distribution
pdfunc(model = "Cauchy", xmin = -4, xmax = 4) 

# Pearson distribution
pdfunc(model = "Pearson", shape = 2, xmin = -4, xmax = 4) 

# exponential distribution
pdfunc(model = "exp", xmin = 0, xmax = 8) 

pdfunc(model = "exp", xmin = -4, xmax = 4)

# Chi-square distribution
pdfunc(model = "Chi2", df = 3, xmin = 0, xmax = 8) 

# double exponential distribution
pdfunc(model = "dexp", xmin = -4, xmax = 2) 

# uniform distribution
pdfunc(model = "unif", xmin = 0, xmax = 1) 
}

\keyword{ts}
