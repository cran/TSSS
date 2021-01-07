\name{klinfo}
\alias{klinfo}

\title{Kullback-Leibler Information}
\description{
  Computes Kullback-Leibler information.}
\usage{
klinfo(distg = 1, paramg = c(0, 1), distf = 1, paramf, xmax = 10)
}
\arguments{
  \item{distg}{function for the true density (1 or 2).
    \tabular{rl}{
      1 : \tab Gaussian (normal) distribution\cr
      \tab paramg(1): mean\cr
      \tab paramg(2): variance\cr
      2 : \tab Cauchy distribution\cr
      \tab paramg(1): \eqn{\mu} (location parameter)\cr
      \tab paramg(2): \eqn{\tau^2} (dispersion parameter)}
  }
  \item{paramg}{parameter vector of true density.}
  \item{distf}{function for the model density (1 or 2).
    \tabular{rl}{
      1 : \tab Gaussian (normal) distribution\cr
      \tab paramf(1): mean\cr
      \tab paramf(2): variance\cr
      2 : \tab Cauchy distribution\cr
      \tab paramf(1): \eqn{\mu} (location parameter)\cr
      \tab paramf(2): \eqn{\tau^2} (dispersion parameter)}
  }
  \item{paramf}{parameter vector of the model density.}
  \item{xmax}{upper limit of integration. lower limit xmin = -\code{xmax}.}
}
\value{
  \item{nint}{number of function evaluation.}
  \item{dx}{delta.}
  \item{KLI}{Kullback-Leibler information, \eqn{I(g;f)}.}
  \item{gint}{integration of \eqn{g(y)} over [-\code{xmax}, \code{xmax}].}
}
\references{
  Kitagawa, G. (2020)
  \emph{Introduction to Time Series Modeling with Applications in R}.
  Chapman & Hall/CRC.
}
\examples{
# g:Gauss, f:Gauss
klinfo(distg = 1, paramg = c(0, 1), distf = 1, paramf = c(0.1, 1.5), xmax = 8)

# g:Gauss, f:Cauchy
klinfo(distg = 1, paramg = c(0, 1), distf = 2, paramf = c(0, 1), xmax = 8)
}

\keyword{ts}
