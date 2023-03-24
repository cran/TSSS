# PROGRAM 7.1 ARFIT
arfit <- function(y, lag = NULL, method = 1, plot = TRUE, ...)
{
  n <- length(y)             # data length
  if (is.null(lag)) {
    lag <- as.integer(2 * sqrt(n))    # highest order of AR model
  } else if (lag < 1) {
    stop("lag is not a positive integer.\n")
  }

# method                     # estimation procedure
  nf <- 200                  # number of frequencies for computing spectrum
# mj2                        # Adjustable dimension in Fortran
  if (method == 2) {
    mj2 <- n
    if (n > 50000)
      mj2 <- 2 * lag
  } else {
    mj2 <- 1
  }

  yy <- y
  ier <- 0
  if (anyNA(yy) == TRUE) {
    if (method == 1) {
      outmax <- 1.0e+30
      yy[is.na(yy)] <- outmax + 1
    } else {
      stop("Cannot handle the data with missing values.\n")
    }
  }

  z <- .Fortran(C_arfit,
                as.double(yy),
                as.integer(n),
                as.integer(lag),
                as.integer(nf),
                as.integer(mj2),
                as.integer(method),
			    sig = double(lag+1),
			    aic = double(lag+1),
		        mar = integer(1),
			    a = double(lag*lag),
			    par = double(lag),
			    sp = double(nf+1))

  mmin <- z$mar
  a <- array(z$a, dim = c(lag, lag))
  arcoef <- list()
  for (i in 1:lag)
    arcoef[[i]] <- a[1:i, i]

  arfit.out <- list(sigma2 = z$sig, maice.order = mmin, aic = z$aic,
                    arcoef = arcoef, parcor = z$par, spec = z$sp,
                    tsname = deparse(substitute(y)))

  class(arfit.out) <- "arfit"

  if (plot) {
    plot.arfit(arfit.out, ...)
    invisible(arfit.out)
  } else arfit.out

}

plot.arfit <- function(x, ...)
{
  lag <- length(x$aic) - 1
  morder <- x$maice.order
  aicmin <- x$aic[morder + 1]

  old.par <- par(no.readonly = TRUE)
  par(mfrow=c(2, 2), xaxs = "i", yaxs = "i")

  ylim <- c(-1, 1)
  xmax <- lag
  mtitle <- c(paste(x$tsname), "Parcor")
  plot(x$parcor, type = "l", xlim = c(0, xmax), ylim = ylim,
       main = mtitle, xlab = "", ylab = "", ...)
  par(new = TRUE)
  plot(x$parcor, type = "h", xlim = c(0, xmax), ylim = ylim, xlab = "lag",
       ylab = "", ...)

  daic <- x$aic - aicmin
  xx <- c(0:lag)
  ymax <- min(max(daic), 50)
  mtitle <- paste("\nminimum AIC = ", format(round(aicmin, 2), nsmall = 2),
                  "(at order", morder, ")")
  ylabel <- paste("aic - aicmin (Truncated at ", ymax, ")")
  if (ymax < quantile(daic, probs = 0.1))
    ymax <- as.integer(quantile(daic, probs = 0.1))
  plot(xx, daic, type = "l", ylim = c(0, ymax), main = mtitle, xlab = "Lag",
       ylab = ylabel, ...)
  par(new = TRUE)
  plot(xx, daic, type = "h", ylim = c(0, ymax), xlab = "", ylab = "", ...)

  nf1 <- length(x$spec)
  xx <- rep(0, nf1)
  for (i in 1:nf1)
    xx[i] <- (i - 1) / (2 * (nf1-1))
  ylim <- c(floor(min(x$spec)), ceiling(max(x$spec)))
  plot(xx, x$spec, type = "l", ylim = ylim, main = "\nPower spectrum",
       xlab = "f", ylab = "log p(f)", ...)

  par(old.par)
}
