# PROGRAM 7.1 ARFIT
arfit <- function(y, lag = NULL, method = 1, plot = TRUE, ...)
{
  n <- length(y)             # data length
  if (is.null(lag))
    lag <- as.integer(2 * sqrt(n))    # highest order of AR model
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

  z <- .Call("ArfitC",
             as.double(y),
             as.integer(n),
             as.integer(lag),
             as.integer(nf),
			 as.integer(mj2),
             as.integer(method))

  sig2 <- z[[1L]]
  aic <- z[[2L]]
  mmin <- z[[3]]
  a <- array(z[[4L]], dim = c(lag, lag))
  arcoef <- a[1:mmin, mmin]
  parcor <- z[[5L]]
  if (method == 2)
    parcor <- parcor[1:mmin]
  spec <- z[[6L]]

  arfit.out <- list(sigma2 = sig2, maice.order = mmin, aic = aic,
                    arcoef = arcoef, parcor = parcor, spec = spec,
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
  plot(x$parcor, type = "l", xlim = c(0, xmax), ylim = ylim,
       main = paste(x$tsname), xlab = "", ylab = "", ...)
  par(new = TRUE)
  plot(x$parcor, type = "h", xlim = c(0, xmax), ylim = ylim, xlab = "lag",
       ylab = "parcor", ...)

  daic <- x$aic - aicmin
  xx <- c(0:lag)
  ymax <- min(max(daic), 50)
  if (ymax < quantile(daic, probs = 0.1))
    ymax <- as.integer(quantile(daic, probs = 0.1))
  plot(xx, daic, type = "l", ylim = c(0, ymax),
     main = paste("aic(", morder, ") = ", format(round(aicmin, 2), nsmall = 2)),
     xlab = "Lag", ylab = paste("aic - aicmin (Truncated at ", ymax, ")"), ...)
  par(new = TRUE)
  plot(xx, daic, type = "h", ylim = c(0, ymax), xlab = "", ylab = "", ...)

  nf1 <- length(x$spec)
  xx <- rep(0, nf1)
  for (i in 1:nf1)
    xx[i] <- (i - 1) / (2 * (nf1-1))
  ylim <- c(floor(min(x$spec)), ceiling(max(x$spec)))
  plot(xx, x$spec, type = "l", ylim = ylim, main = "Power spectrum in log scale",
       xlab = "f", ylab = "log p(f)", ...)

  par(old.par)
}
