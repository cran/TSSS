# PROGRAM 11.1
polreg <- function (y, order, plot = TRUE, ...) 
{
  n <- length(y)      # data length
  k <- order          # order of polinomial regression  
  k1 <- k + 1

  z <- .Fortran(C_polreg,
                as.double(y),
                as.integer(n),
                as.integer(k),
                a = double(k * k),
                sig2 = double(k1),
                aic = double(k1),
                imin = integer(1),
                data = double(n))

  a <- array(z$a, dim = c(k, k))
  coef <- list()
  for (i in 1:k)
    coef[[i]] <- a[1:i, i]
  aic <- z$aic
  imin <- z$imin
  daic <- aic - aic[imin + 1]

  polreg.out <- list(order.maice = imin, sigma2 = z$sig2, aic = aic, daic = daic,
                     coef = coef, trend = z$data)
  class(polreg.out) <- "polreg"

  if (plot) {
    rdata <- deparse(substitute(y))
    eval(parse(text=paste(rdata, "<- y")))
    eval(parse(text=paste("plot.polreg(polreg.out,", rdata, ", ...)")))
    invisible(polreg.out)
  } else polreg.out
}

plot.polreg <- function(x, rdata = NULL, ...)
{
  ts.atr <- tsp(rdata)
  trend <- x$trend
  if (is.null(ts.atr) == FALSE)
    trend <- ts(trend, start = ts.atr[1], frequency = ts.atr[3])

  old.par <- par(no.readonly = TRUE)

  morder <- x$order.maice
  aicmin <- format(round(x$aic[morder+1], 2), nsmall=2)
  ylim <- range(trend, rdata, na.rm = TRUE)

  if (is.null(rdata) == TRUE) {
    mtitle <- paste("Trend component\n minimum aic =",
                    aicmin, " at order", morder)
  } else {
    tsname <- deparse(substitute(rdata))
    mtitle <- paste(tsname, "and trend component\n minimum aic =", aicmin,
                    " at order", morder)
    plot(rdata, type = "l", ylim = ylim, xlab = "", ylab = "", main = "",
         xaxs = "i", xaxt = "n", ...)
    par(new = TRUE)
  }

  plot(trend, type = "l", ylim = ylim, col = 2, xlab = "", main = mtitle,
       ylab = "", xaxs = "i", cex.main = 1.0, ...) 
  par(old.par)
}