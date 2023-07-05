# PROGRAM 5.5 LSQR
lsqr <- function(y, lag = NULL, period = 365, plot = TRUE, ...)
{
  n <- length(y)                 # length of data
  if (is.null(lag))
    lag <- as.integer(sqrt(n))   # number of sine and cosine terms
  if (as.integer(period/2) < lag) {
    lag <- as.integer(period/2)
    warning(gettextf("lag is corrected to %d. \n
            Maximum number of lag should be less than or equal to floor(period/2).", lag), domain = NA)
  }
  k <- 2 * lag + 1

  mj <- n
  if (n > 50000) mj <- 2 * k

  z <- .Fortran(C_lsqr,
                as.double(y),
                as.integer(n),
                as.integer(k),
                as.integer(period),
                as.integer(mj),
                aic = double(k + 1),
                sig2 = double(k + 1),
                imin = integer(1),
                reg = double(k * k),
                data = double(n))

  a <- array(z$reg, dim = c(k, k))
  reg <- list()
  for (i in 1:k)
    reg[[i]] <- a[1:i, i]

  lsqr.out <- list(aic = z$aic, sigma2 = z$sig2, maice.order = z$imin,
                   regress = reg, tripoly = z$data)
  class(lsqr.out) <- "lsqr"

  if (plot) {
    rdata <- deparse(substitute(y))
    eval(parse(text=paste(rdata, "<- y")))
    eval(parse(text=paste("plot.lsqr(lsqr.out,", rdata, ", ...)")))
  }

  lsqr.out
}

print.lsqr <- function(x, ...)
{
  aic <- x$aic
  nl <- length(aic)
  sig2 <- x$sigma2
  morder <- x$maice.order
  aicmin <- aic[morder+1] 

  message("\n Order     Sigma2            AIC         AIC-AICMIN")
  for (i in 1:nl) 
    message(gettextf("%5i  %13.6e   %12.3f   %12.3f", i-1, sig2[i], aic[i],
            aic[i] - aicmin), domain = NA)

  message(gettextf("\n Minimum AIC = %12.3f   attained at m = %5i", aicmin,
          morder), domain = NA)
  message(gettextf(" Sigma2 = %13.6e\n", sig2[morder+1]), domain = NA)
}

plot.lsqr <- function(x, rdata = NULL, ...)
{
  ts.atr <- tsp(rdata)
  data <- x$tripoly
  if (is.null(ts.atr) == FALSE)
    data <- ts(data, start = ts.atr[1], frequency = ts.atr[3])
  imin <- x$maice.order

  if (is.null(rdata) == FALSE) {
    ylim <- c(floor(min(data, rdata)), ceiling(max(data, rdata)))
    tsname <- deparse(substitute(rdata))
    mtitle <-
      paste(tsname, "\nand regression curve of the model with order", imin)
    plot(rdata, type = "l", ylim = ylim, xlab = "", ylab = "", xaxs = "i", ...)
    par(new = TRUE) 
  } else {
    ylim <- c(floor(min(data)), ceiling(max(data)))
    mtitle <-
      paste("Rregression curve of the model with order",imin)
  }
  plot(data, type = "l", ylim = ylim, col=2, xlab = "", ylab = "",
       main = mtitle, xaxs = "i", cex.main = 1.0, lwd = 2, ...)

}
