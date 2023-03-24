\name{HAKUSAN}
\docType{data}
\alias{HAKUSAN}
\title{Ship's Navigation Data}
\description{
  A multivariate time series of a ship's yaw rate, rolling, pitching and rudder
  angles which were recorded every second while navigating across the Pacific
  Ocean.
}
\usage{
data(HAKUSAN)
}
\format{
  A data frame with 1000 observations on the following 4 variables.
  \tabular{rll}{
  [, 1] \tab YawRate \tab yaw rate\cr
  [, 2] \tab Rolling \tab rolling\cr
  [, 3] \tab Pitching\tab pitching\cr
  [, 4] \tab Rudder  \tab rudder angle
  }
}
\source{
  The data were offered by Prof. K. Ohtsu of Tokyo University of Marine Science
  and Technology.
}
\examples{
data(HAKUSAN)
HAKUSAN234 <- HAKUSAN[, c(2,3,4)]

## put histograms on the diagonal
panel.hist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
#    par(usr = c(usr[1:2], 0, 1.4))
    par(usr = c(usr[1:2], 0, 1.3))
    nB <- 20; nB1 <- nB + 1
    xmin <- min(x)
    xmax <- max(x)
    w <- (xmax - xmin) / nB
    breaks <- xmin
    b <- xmin
    for (i in 1:nB) {
      b <- b + w
      breaks <- c(breaks, b)
    }
    h <- hist(x, breaks = breaks, plot = FALSE)
    y <- h$counts; y <- y / max(y)
    rect(breaks[1:nB], 0, breaks[2:nB1], y, ...)

}

par(xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n")
pairs(HAKUSAN234, diag.panel = panel.hist, pch = 20, cex.labels = 1.5, 
      label.pos = 0.9, lower.panel = NULL)
}
\keyword{datasets}
