# PROGRAM 11.2
trend <- function (y, trend.order = 1, tau2.ini = NULL, delta, plot = TRUE, ...)     
{
  n <- length(y)      # data length
  m <- trend.order

# tau2.ini                 # initial variance of systen noise
# delta                  # delta for computing variance of system noise
  iopt <- 1               # search method
  if (is.null(tau2.ini)) {
    iopt  <- 0
    tau2.ini <- 0
    delta <- 0
  }

  z <- .Call("TrendC",
             as.double(y),
             as.integer(n),
             as.integer(m),
             as.integer(iopt),
             as.double(tau2.ini),
             as.double(delta))

  xss <- array(z[[5L]], dim = c(m, n))
  trend <- xss[1, 1:n]
  res <- z[[6L]]

  trend.out <- list(trend = trend, residual = res, tau2 = z[[1L]],
                    sigma2 = z[[2L]], llkhood = z[[3L]], aic = z[[4L]])
  class(trend.out) <- "trend"

  if(plot) {
    rdata <- deparse(substitute(y))
    eval(parse(text=paste(rdata, "<- y")))
    eval(parse(text=paste("plot.trend(trend.out,", rdata, ", ...)")))
  }

  return(trend.out)
}

print.trend <- function(x, ...)
{
  message(gettextf("\n tau2 \t\t%12.5e", x$tau2), domain = NA)
  message(gettextf(" sigma2 \t%12.5e", x$sigma2), domain = NA)
  message(gettextf(" log-likelihood\t%12.3f", x$llkhood), domain = NA)
  message(gettextf(" aic\t\t%12.3f\n", x$aic), domain = NA)
}

plot.trend <- function(x, rdata = NULL, ...)
{
    ts.atr <- tsp(rdata)
    trend <- x$trend
    res <- x$res
    if (is.null(ts.atr) == FALSE) {
      trend <- ts(trend, start = ts.atr[1], frequency = ts.atr[3])
      res <- ts(res, start = ts.atr[1], frequency = ts.atr[3])
    }

    old.par <- par(no.readonly = TRUE)
    par(mfcol = c(2, 1), xaxs = "i")

    ylim <- range(x$trend, rdata, na.rm = TRUE)
    if (is.null(rdata) == TRUE) {
      mtitle <- paste("Trend component")
    } else {
      tsname <- deparse(substitute(rdata))
      mtitle <- paste(tsname, "and trend component")
      plot(rdata, type = "l", ylim = ylim, xlab = "", ylab = "", main = "", ...) 
      par(new = TRUE)
    }
    plot(trend, type = "l", ylim = ylim, col=2, xlab = "", ylab = "",
         main = mtitle, ...) 

    plot(res, type = "h", xlab = "",  ylab = "",  main = "Residuals", ...)
    abline(h = 0)
    par(old.par)
}