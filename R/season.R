# PROGRAM 12.1
season <- function(y, trend.order = 1, seasonal.order = 1, ar.order = 0,
                   trade = FALSE, period = 12, tau2.ini = NULL,
                   filter = c(1,length(y)), predict = length(y),
                   arcoef.ini = NULL, log = FALSE,
                   minmax = c(-1.0e+30, 1.0e+30), plot = TRUE, ...) 
{
  if (trend.order < 0)
    stop("trend.order must be greater than 0 or equal to 0.")

#  if ((trade != 0) && (trade != 6 ))
#    stop("'trade' is 0 or 6" )

#  if (trade == 6 )
#    if (is.null(year) )
#      stop("specify starting year of the data" )

  if (is.null(tsp(y)) == TRUE) {
	year <- 1
	month <- 1
  } else if (is.null(tsp(y)) == FALSE) {
	year <- start(y)[1]
	month <- start(y)[2]
    period <- tsp(y)[3]
  }


  n <- length(y)      # data length
  m1 <- trend.order
  m2 <- seasonal.order
  m3 <- ar.order
#  m4 <- trade
  m4 <- 0
  if (trade == TRUE) {
    if(period == 4 || period ==12) {
      m4 <- 6
    } else {
      warning("'trade = TRUE' is available only if period is 4 or 12.")
    }
  }
  jyear <- year
#  if (is.null(jyear) ) jyear <- 0
  jmonth <- month
  logt <- 0
  if (log == TRUE) logt <- 1
  iopt <- 1
  outmin <- minmax[1]
  outmax <- minmax[2]
  ns <- filter[1]
  nfe <- filter[2]
  npe <- predict
  nmax <- max(n, npe)

  np <- 0
  mj <- 0
  if (m1 > 0) {
    np <- 1
    mj <- m1
  }
  if (m2 > 0) {
    np <- np+1
    mj <- mj + m2*(period-1)
  }
  if (m3 > 0) {
    np <- np+1
    mj <- mj + m3
  }
  if (m4 > 0)
    mj <- mj + 6

  tau2 <- rep(0, 4)
  ll <- 0
  if (is.null(tau2.ini) ) {
    if (m1 > 0 ) {
      ll <- 1
      tau2[1] <- 0.05
    }
    if (m2 > 0 ) {
      ll <- ll+1
      tau2[ll] <- 0.11e-7
    }
    if (m3 > 0 ) {
      ll <- ll+1
      tau2[ll] <- 0.9
    }
    if (m4 > 0 ) {
      ll <- ll+1
      tau2[ll] <- 0.0
    }
  } else {
    nc <- length(tau2.ini)
    tau2[1:nc] <- tau2.ini[1:nc]
    for (i in 1:nc)
      if (tau2[i] == 1)
        stop("variance of the system noise tau2 is not equal to 1")
  }

  arcoef <- 0
  if (m3 > 0) {
    if (is.null(arcoef.ini)) {
      par <- rep(0,m3)
      a <- rep(0,m3)
      arcoef <- rep(0,m3)
      for (i in 1:m3)
        par[i] <- -(-0.8)**i
      for (i in 1:m3 ) {
        arcoef[i] <- par[i]
        a[i] <- par[i]
        if (i > 1 ) {
          im1 <- i-1
          for (j in 1:im1)
            arcoef[j] <- a[j] - par[i]*a[i-j] 
          if (i < m3 ) for (j in 1: im1 ) a[j] <- arcoef[j]
        }
      }
    } else {
      if (length(arcoef.ini) < m3 )
        stop(gettextf( "%d specify initial AR coefficients", m3 ), domain = NA)
      arcoef <- arcoef.ini
    }
  }

  z <- .Call("SeasonC",
             as.double(y),
             as.integer(n),
             as.integer(m1),
             as.integer(m2),
             as.integer(m3),
             as.integer(m4),
             as.integer(period),
             as.integer(jyear),
             as.integer(jmonth),
             as.double(tau2),
             as.integer(ns),
             as.integer(nfe),
             as.integer(npe),
             as.double(arcoef),
             as.integer(logt),
             as.integer(iopt),
             as.double(outmin),
             as.double(outmax),
             as.integer(nmax),
             as.integer(mj))

  ier1 <- z[[7L]]
  if (ier1 != 0)
    stop("Log-transformation cannot be applied to zeros and nagative numbers")

  ier2 <- z[[8L]]
  if (ier2 == 1)
    stop(" Matrix with zero row in decompose" )
  if (ier2 == 2)
    stop(" Singular matrix in decompose.zero divide in solve" )
  if (ier2 == 3)
    stop(" Convergence in impruv.matrix is nearly singular" )

  trend <- NULL
  seasonal <- NULL
  ar <- NULL
  if (m3 == 0) arcoef <- NULL
  deff <- NULL

  sig2 <- z[[2L]]
  m12 <- m1 + m2*(period-1)
  xss <- array(z[[4L]], dim = c(mj, nmax))
  if (m1 > 0)
    trend <- xss[1,1:nmax]
  if (m2 > 0)
    seasonal <- xss[m1+1,1:nmax]
  if (m3 > 0)
    ar <- xss[m12+1, 1:nmax]
  if (m4 > 0)
    deff <- z[[6L]]
  vss <- array(z[[5L]], dim = c(mj,mj,nmax))
  vss <- vss[1,1,1:nmax]

  yy <- y
  if (log == TRUE)
    yy <- log10(y)

  noise <- yy
  if (m1 > 0)
    noise <- noise - trend[1:n]
  if (m2 > 0)
    noise <- noise - seasonal[1:n]
  if (m3 > 0)
    noise <- noise - ar[1:n]
  if (m4 > 0)
    noise <- noise - deff[1:n]
  noise.list <- list(val = noise, range = c(ns, nfe))

  season.out <- list(tau2=tau2[1:np], sigma2=sig2, llkhood=z[[1L]], aic=z[[3L]],
                     trend=trend, seasonal=seasonal, arcoef=arcoef, ar=ar,
                     day.effect=deff, noise=noise.list, cov = vss)
  class(season.out) <- "season"

  if (plot == TRUE) {
    rdata <- deparse(substitute(y))
    eval(parse(text=paste(rdata, "<- y")))
    if (log == TRUE) {
      eval(parse(text=paste("plot.season(season.out, log10(", rdata, "), ...)")))
    } else {
      eval(parse(text=paste("plot.season(season.out,", rdata, ", ...)")))
    }
  }

  return(season.out)
}

print.season <- function(x, ...)
{
  n <- length(x$tau2)
  message("\n tau2\t", appendLF=FALSE)
  for (i in 1:n)
    message(gettextf("\t%12.5e", x$tau2[i]), appendLF=FALSE, domain = NA)
  message(gettextf("\n sigma2\t\t%12.5e", x$sigma2), domain = NA)
  message(gettextf(" log-likelihood\t%12.3f", x$llkhood), domain = NA)
  message(gettextf(" aic\t\t%12.3f\n", x$aic), domain = NA)
}

plot.season <- function(x, rdata = NULL, ...)
{
  ts.atr <- tsp(rdata)
  n <- length(rdata)
  trend <- x$trend
  seasonal <- x$seasonal
  ar <- x$ar
  day.effect <- x$day.effect
  noise <- x$noise$val
  if (is.null(ts.atr) == FALSE)
    noise <- ts(noise, start = ts.atr[1], frequency = ts.atr[3])

  y <- rdata
  n <- length(noise)
  range <- x$noise$range
  if (is.null(range) == TRUE) {
    ns <- 1
    nfe <- n
  } else {
    ns <- range[1]
    nfe <- range[2]
  }
  npe <- max(length(trend), length(seasonal), length(ar), length(day.effect))

  np <- length(x$tau2)
  ng <- np + 1
  if (is.null(day.effect) == FALSE)
    ng <- ng + 1
  nc <- 1
  nr <- ng
  if (nfe < npe)  ng <- ng + 1
  if (ng > 3) {
    nc <- 2
    nr <- as.integer((ng + 1) / 2)
  }

  old.par <- par(no.readonly = TRUE)
  par(mfrow = c(nr, nc), xaxs = "i")

  nmax <- max(n, npe)
  xlim <- c(1, nmax)
  if (is.null(ts.atr) == FALSE)
    xlim <- xlim + ts.atr[1] - 1
  ylim1 <- NULL

### original and trend
  if (is.null(trend) == FALSE) {
    trendx <- c(ns:nmax)
    if (is.null(ts.atr) == FALSE) {
      trend <- ts(trend, start = ts.atr[1], frequency = ts.atr[3])
      trendx <- trendx + ts.atr[1] - 1
    }
    s2 <- x$sigma2
    t1 <- trend[1:nmax]
    t2 <- trend[1:nmax]
    for (i in 1:nmax) {
      t1[i] <- t1[i] - sqrt(x$cov[i] * s2)
      t2[i] <- t2[i] + sqrt(x$cov[i] * s2)
    }
    ylim1 <- range(t1[ns:nmax], t2[ns:nmax], y)
    if (is.null(y) == TRUE) {
      mtitle <- paste("Trend component")
    } else {
      tsname <- deparse(substitute(rdata))
      mtitle <- paste(tsname, "and trend component")
      plot(y, type = "l", ylim = ylim1, xlab = "", ylab = "", main = mtitle, ...)
      par(new = TRUE)
    }
    plot(trendx, trend[ns:nmax], type = "l", col = 2, xlim = xlim, ylim = ylim1, xlab = "",
         ylab = "", main = mtitle, xaxt = "n", ...) 
    par(new = TRUE)
    plot(trendx, t1[ns:nmax], type = "l", lty = 3, col = 4, xlim = xlim, ylim = ylim1, xlab = "",
         ylab = "", main = "", xaxt = "n", ...)
    par(new = TRUE)
    plot(trendx, t2[ns:nmax], type = "l", lty = 3, col = 4, xlim = xlim, ylim = ylim1, xlab = "",
         ylab = "", main = "", xaxt = "n", ...)
  }

  ylim <- range(seasonal, ar, day.effect, na.rm = TRUE)
  ymax <- max(abs(ylim[1]), abs(ylim[2]))
  ylim <- c(-ymax, ymax)

### seasonal
  if (is.null(seasonal) == FALSE) {
    if (is.null(ts.atr) == FALSE)
      seasonal <- ts(seasonal, start = ts.atr[1], frequency = ts.atr[3])
    plot(seasonal, type = "h", ylim = ylim, xlab = "", ylab = "",
         main = "Seasonal component", ...) 
  }

### AR
  if (is.null(x$ar) == FALSE) {
    ar <- x$ar
    if (is.null(ts.atr) == FALSE)
      ar <- ts(ar, start = ts.atr[1], frequency = ts.atr[3])
    plot(ar, type = "h", xlab = "", ylim = ylim, ylab = "",
         main = "AR component", ...)
  }

### day.effect
  if (is.null(x$day.effect) == FALSE) {
    day.effect <- x$day.effect
    if (is.null(ts.atr) == FALSE)
      day.effect <- ts(day.effect, start = ts.atr[1], frequency = ts.atr[3])
    plot(day.effect, type = "h", xlab = "", ylim = ylim, ylab = "",
         main = "Trading day effect", ...) 
  }

### noise
  noisex <- c(ns:nfe)
  if (is.null(ts.atr) == FALSE) {
    noise <- ts(noise, start = ts.atr[1], frequency = ts.atr[3])
    noisex <- noisex + ts.atr[1] - 1
  }
  plot(noise, type = "n", xlab = "", ylim = ylim, ylab = "", main = "noise", ...)
  par(new = TRUE)
  plot(noisex, noise[ns:nfe], type = "h", xlab = "", xlim = xlim, ylim = ylim,
       ylab = "", xaxt = "n", ...) 

  if (nfe < npe) {
    nps <- nfe + 1
    yy <- rep(0, npe)
    for (i in nps:npe) {
      if (is.null(seasonal) == FALSE)
        yy[i] <- yy[i] + seasonal[i]
      if (is.null(ar) == FALSE)
        yy[i] <- yy[i] + ar[i]
      if (is.null(day.effect) == FALSE)
        yy[i] <- yy[i] + day.effect[i]
    }

    xlim <- c(0, nmax)
    xx <- c(nps:npe)
    npsv <- nps
    if (is.null(ts.atr) == FALSE) {
      xlim <- xlim + ts.atr[1] - 1
      xx <- xx + ts.atr[1] - 1
      npsv <- npsv + ts.atr[1] - 1
    }

    if (is.null(trend) == FALSE) {
      for (i in nps:npe)
        yy[i] <- yy[i] + trend[i]
      t1 <- yy
      t2 <- yy
      for (i in nps:npe) {
        t1[i] <- t1[i] - sqrt(x$cov[i] * s2)
        t2[i] <- t2[i] + sqrt(x$cov[i] * s2)
      }
      ylim <- range(y, ylim1, t1[nps:npe], t2[nps:npe])
      plot(y, type = "l", ylim = ylim, main = "", xlab = "", ylab = "", ...) 
      par(new = TRUE)
      plot(xx, yy[nps:npe], type = "l", col = 2, xlim = xlim, ylim = ylim, xlab = "",
           ylab = "", main = "predicted values", xaxt = "n", ...)
      par(new = TRUE)
      plot(xx, t1[nps:npe], type = "l", col = 4, xlim = xlim, ylim = ylim, xlab = "",
           ylab = "", main = "", xaxt = "n", ...) 
      par(new = TRUE)
      plot(xx, t2[nps:npe], type = "l", col = 4, xlim = xlim, ylim = ylim, xlab = "",
           ylab = "", main = "", xaxt = "n", ...) 
      abline(v = npsv, lty = 2)

    } else {
      plot(y, type = "l", ylim = ylim, xlab = "", ylab = "",
           main = "", ...) 
      par(new = TRUE)
      plot(xx, yy, type = "l", col = 4, ylim = ylim, xlab = "",
           ylab = "", main = "", ...) 
      abline(v = npsv, lty = 2)
    }
  }
  par(old.par)
}
