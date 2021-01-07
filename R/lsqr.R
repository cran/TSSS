# PROGRAM 5.5 LSQR
lsqr <- function(y, lag = NULL, plot = TRUE, ...)
{
  n <- length(y)                 # length of data
  if (is.null(lag))
    lag <- as.integer(sqrt(n))   # number of sine and cosine terms
  k <- 2 * lag + 1

  mj <- n
  if (n > 50000) mj <- 2 * k

  z <- .Call("LsqrC",
             as.double(y),
             as.integer(n),
             as.integer(k),
             as.integer(mj))

  aic <- z[[1L]]
  sig2 <- z[[2L]]
  imin <- z[[3L]]
  a <- array(z[[4L]], dim = c(k, k))
  data <- z[[5L]]
  reg <- list()
  for (i in 1:k)
    reg[[i]] <- a[1:i, i]

  lsqr.out <- list(aic = aic, sigma2 = sig2, maice.order = imin, regress = reg,
                   tripoly = data)
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
       main = mtitle, xaxs = "i", cex.main = 1.0,...)

}

