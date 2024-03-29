# PROGRAM 2.1
unicor <- function (y, lag = NULL, minmax = c(-1.0e+30, 1.0e+30), plot = TRUE, ...)
{
  n <- length(y)
  if (is.null(lag))
    lag <- as.integer(2 * sqrt(n))    # maximum lag
  lag1 <- lag + 1
  outmin <- minmax[1]
  outmax <- minmax[2]

  z <- .Fortran(C_unicorf,
                as.double(y),
                as.integer(n),
                as.integer(lag),
                as.double(outmin),
                as.double(outmax),
                cov = double(lag1 * 4),
                mean = double(1))

  cov <- matrix(z$cov, lag1, 4)
  acov <- cov[, 1]
  acor <- cov[, 2]
  cerr <- cov[, 3]
  rerr <- cov[, 4]

  unicor.out <- list(acov = acov, acor = acor, acov.err = cerr, acor.err = rerr,
                     mean = z$mean, tsname = deparse(substitute(y)))
  class(unicor.out) <- "unicor"

  if (plot) {
    plot.unicor(unicor.out, ...)
    invisible(unicor.out)
  } else unicor.out
}

plot.unicor <- function(x, ...) {
  acor <- x$acor
  lag <- length(acor)[1] - 1

  old.par <- par(no.readonly = TRUE)
  mtitle <- paste(x$tsname)
  plot(c(0:lag), acor, type = "h", ylim = c(-1, 1), main = mtitle,
       ylab = "Autocorrelation", xlab = "Lag", xaxs = "i", ...)
  par(new = TRUE)
  plot(c(0:lag), acor, type = "l", ylim = c(-1, 1), ylab = "", xlab = "",
       xaxs = "i", ...)
  par(new = TRUE)
  abline(h = 0)
  par(old.par)
}
