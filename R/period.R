# PROGRAM 3.1
period <-
  function(y, window = 1, lag = NULL, minmax = c(-1.0e+30, 1.0e+30),
           plot = TRUE, ...) {

  n <- length(y)

  if (window == 0) {
    if (is.null(lag))
      lag <- n - 1
    np <- (lag + 1) /2
  } else if (window == 1 || window == 2) {
    if (is.null(lag))
      lag <- as.integer(2 * sqrt(n))
    np <- lag
  } else {
    stop("'window' is 0, 1 or 2")
  }

  outmin <- minmax[1]
  outmax <- minmax[2]

  z <- .Call("PeriodC",
             as.double(y),
             as.integer(n),
             as.integer(np),
             as.integer(window),
             as.integer(lag),
             as.double(outmin),
             as.double(outmax))

  np1 <- np + 1
  pe <- z[[1]][1:np1]     # periodogram (raw spectrum)
  spe <- z[[2]][1:np1]    # logarithm of smoothed periodogram
  ifg <- z[[4L]]
  if (ifg == 0) log <- "TRUE"
  if (ifg != 0) log <- "FALSE"

  period.out <- list(period = pe, smoothed.period = spe, window = window,
                     log.scale = log, tsname = deparse(substitute(y)))
  class(period.out) <- "spg"

  if (plot) {
    plot.spg(period.out, ...)
    invisible(period.out)
  } else period.out
}


# PROGRAM 3.2
fftper <- 
  function(y, window = 1, plot = TRUE, ...) {

  n <- length(y)
  if ((window < 0) || (window > 2 ))
    stop("'window' is 0, 1 or 2")

  z <- .Call("FftperC",
             as.double(y),
             as.integer(n),
             as.integer(window))

  np <- z[[3]]
  np1 <- np + 1
  pe <- z[[1]][1:np1]     # periodogram (raw spectrum)
  spe <- z[[2]][1:np1]    # logarithm of smoothed periodogram
  ifg <- z[[4L]]
  if (ifg == 0) log <- "TRUE"
  if (ifg != 0) log <- "FALSE"

  fftper.out <- list(period = pe, smoothed.period = spe, window = window,
                     log.scale = log, tsname = deparse(substitute(y)))
  class(fftper.out) <- "spg"

  if (plot) {
    plot.spg(fftper.out, ...)
    invisible(fftper.out)
  } else fftper.out
}

plot.spg <- function(x, type = "vl", ...)
{
  spe <- x$smoothed.period
  np1 <- length(spe)
  np <- np1 - 1
  spe <- spe[2:np1]
  xx <- rep(0, np)
  for (i in 1:np)
    xx[i] <- i / (np * 2)
  ymin <- floor(min(spe))
  ylim <- c(ymin, ceiling(max(spe)))

  mtitle <- paste(x$tsname)
  if (x$log.scale == "TRUE") {
    if (x$window == 0)
      ylabel <- "log-periodgram"
    if (x$window != 0)
      ylabel <- "smoothed log-periodgram"
  } else {
    if (x$window == 0)
      ylabel <- "periodgram"
    if (x$window != 0)
      ylabel <- "smoothed periodgram"
  }

  if (type == "l") {
    plot(xx, spe, type = "l", xlim = c(0, 0.5), main = mtitle, ylim = ylim,
         xlab = "f", ylab = ylabel, ...)
  } else {
    plot(xx, spe, type = "n", xlim = c(0, 0.5), main = mtitle, ylim = ylim,
         xlab = "f", ylab = ylabel, ...)
    for (i in 1:np)
      segments(xx[i], ymin, xx[i], spe[i], lwd = 0.7, ...)
  }
  
}

