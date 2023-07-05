# PROGRAM 11.2
trend <- function (y, trend.order = 1, tau2.ini = NULL, delta, plot = TRUE, ...)     
{
  n <- length(y)      # data length
  m <- trend.order
  if (m < 1)
    stop("Trend order must be a positive integer.")


# tau2.ini                 # initial variance of systen noise
# delta                  # delta for computing variance of system noise
  iopt <- 1               # search method
  if (is.null(tau2.ini)) {
    iopt  <- 0
    tau2.ini <- 0
    delta <- 0
  }

  z <- .Fortran(C_trend,
                as.double(y),
                as.integer(n),
                as.integer(m),
                as.integer(iopt),
                as.double(tau2.ini),
                as.double(delta),
                tau2 = double(1),
                sig2 = double(1),
                lkhood = double(1),
                aic = double(1),
                xss = double(n * m),
                vss = double(n * m * m),
                rs = double(n))

  xss <- array(z$xss, dim = c(m, n))
  trend <- xss[1, 1:n]
  vss <- array(z$vss, dim = c(m,m,n))
  vss <- vss[1,1,1:n]
  res <- z$rs

  ts.atr <- tsp(y)
    if (is.null(ts.atr) == FALSE) {
      trend <- ts(trend, start = ts.atr[1], frequency = ts.atr[3])
      res <- ts(res, start = ts.atr[1], frequency = ts.atr[3])
    }

  trend.out <- list(trend = trend, residual = res, tau2 = z$tau2,
                    sigma2 = z$sig2, llkhood = z$lkhood, aic = z$aic, cov = vss)
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

    trend <- x$trend
    res <- x$res
    ts.atr <- tsp(trend)

    if (is.null(ts.atr) == TRUE) {
      xtime <- "time"
    } else {
      freq <- ts.atr[3]
      if (freq == 4 || freq == 12) {
        xtime <- "year"
      } else if (freq == 24) {
        xtime <- "day"
      } else if (freq == 52 || freq == 365/7) {
        xtime <- "week"
      } else if (freq == 365.25/7 || freq == 52.18) {
        xtime <- "week"
      } else {
        xtime <- "time"
      }
    }

    old.par <- par(no.readonly = TRUE)
    par(mfcol = c(2, 1), xaxs = "i")

    ylim <- range(x$trend, rdata, na.rm = TRUE)
    if (is.null(rdata) == TRUE) {
      mtitle <- paste("Trend component")
    } else {
      tsname <- deparse(substitute(rdata))
      mtitle <- paste(tsname, "and trend component")
      plot(rdata, type = "l", ylim = ylim, xlab = "", ylab = "", main = "", xaxt = "n", ...)
      par(new = TRUE)
    }
    plot(trend, type = "l", col = 2, ylim = ylim, xlab = xtime,
         ylab = "", main = mtitle, ...)

    plot(res, type = "h", xlab = xtime,  ylab = "",  main = "Residuals", ...)
    abline(h = 0)

    par(old.par)
}