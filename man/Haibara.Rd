\name{Haibara}
\docType{data}
\alias{Haibara}
\title{Haibara Data}
\description{
A bivariate time series of the groundwater level and the atmospheric pressure
 that were observed at 10-minuite intervals at the Haibara observatory of the
 Tokai region, Japan.
}
\usage{
data(Haibara)
}
\format{
A data frame with 400 observations on the following 2 variables.
  \tabular{rl}{
  [, 1] \tab Groundwater level\cr
  [, 2] \tab Atmospheric pressure
  }
}
\source{
The data were offered by Dr. M. Takahashi and Dr. N. Matsumoto of National
 Institute of Advanced Industrial Science and Technology.
}
\examples{
data(Haibara)

## put histograms on the diagonal
panel.hist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
#    par(usr = c(usr[1:2], 0, 1.4))
    par(usr = c(usr[1:2], 0, 1.3))
    nB <- 15; nB1 <- nB + 1
    xmin <- min(x, na.rm = TRUE)
    xmax <- max(x, na.rm = TRUE)
    w <- (xmax - xmin) / nB
    breaks <- xmin
    b <- xmin
    for (i in 1:nB) {
      b <- b + w
      breaks <- c(breaks, b)
    }
    h <- hist(x, breaks = breaks, plot = FALSE)
    y <- h$counts
    y <- y / max(y)
    rect(breaks[1:nB], 0, breaks[2:nB1], y, ...)
}

par(xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n")
pairs(Haibara, diag.panel = panel.hist, pch = 20, cex.labels = 1.5,
      label.pos = 0.9, lower.panel = NULL)
}
\keyword{datasets}
