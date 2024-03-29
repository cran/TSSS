# PROGRAM 12.1
season <- function(y, trend.order = 1, seasonal.order = 1, ar.order = 0,
                   trade = FALSE, period = NULL, tau2.ini = NULL,
                   filter = c(1,length(y)), predict = length(y),
                   arcoef.ini = NULL, log = FALSE, log.base = "e",
                   minmax = c(-1.0e+30, 1.0e+30), plot = TRUE, ...) 
{
  if (trend.order < 0 || trend.order > 3)
    stop("'trend.order' must be 0, 1, 2 or 3.")

  if (seasonal.order < 0 || seasonal.order > 2)
    stop("'seasonal.order' must be 0, 1 or 2.")

  if (ar.order < 0 || ar.order > 5)
    stop("'ar.order' must be 0, 1, 2, 3, 4 or 5.")

  n <- length(y)      # data length
  m1 <- trend.order
  m2 <- seasonal.order
  m3 <- ar.order
  m4 <- 0
  m123 <- m1 + m2 + m3

  if (m123 == 0)
    stop("Total order is 0.")

  if (m2 !=0 || trade == TRUE)
    if (is.null(tsp(y)) == TRUE && is.null(period) == TRUE)
      stop("Neither 'frequency' and 'period' is not specified.")

  period.in <- period
  if (is.null(tsp(y)) == TRUE) {
	year <- 0
	month <- 1
  } else if (is.null(tsp(y)) == FALSE) {
	year <- start(y)[1]
	month <- start(y)[2]
    period <- tsp(y)[3]
    if (is.null(period.in) == FALSE && period != period.in) {
      warning(gettextf("The tsp attribute is not NULL, the value of 'period' is ignored.
                       \n'period' was set to 'frequency' %d.", period),
                       domain=NA)
    }
  }

  if (m2 != 0) {
    if (is.element(period, c(4, 5, 7, 12, 24)) == FALSE) {
      warning("'period' must be 4, 5, 7, 12 or 24.\n")
      m2 <- 0
    }
  }

  if (trade == TRUE) {
    if(period == 4 || period ==12) {
      m4 <- 6
    } else {
      warning("In the case of trade = TRUE, 'period' must be 4 or 12.\n")
    }
  }

  jyear <- year
  jmonth <- month
  logt <- 0
  if (log == TRUE) {
    if (log.base == "10") {
      logt <- 1
    } else if (log.base == "e") {
      logt <- 2
    } else {
      warning('log.base must be "e" (default) or "10".\n')
    }
  }
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
    if (m1 > 0) {
      ll <- 1
      tau2[1] <- 0.05
    }
    if (m2 > 0) {
      ll <- ll+1
      tau2[ll] <- 0.11e-7
    }
    if (m3 > 0) {
      ll <- ll+1
      tau2[ll] <- 0.9
    }
    if (m4 > 0) {
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

  z <- .Fortran(C_season,
                as.double(y),
                as.integer(n),
                as.integer(m1),
                as.integer(m2),
                as.integer(m3),
                as.integer(m4),
                as.integer(period),
                as.integer(jyear),
                as.integer(jmonth),
                tau2 = as.double(tau2),
                as.integer(ns),
                as.integer(nfe),
                as.integer(npe),
                as.double(arcoef),
                as.integer(logt),
                as.integer(iopt),
                as.double(outmin),
                as.double(outmax),
                as.integer(nmax),
                as.integer(mj),
                ff = double(1),
                over = double(1),
                aic = double(1),
                xss = double(nmax * mj),
                vss = double(nmax * mj * mj),
                deff = double(nmax),
                ier1 = integer(1),
                ier2 = integer(1))

  ier1 <- z$ier1
  if (ier1 == -1) 
    stop("Log-transformation cannot be applied to zeros and nagative numbers")

  ier2 <- z$ier2
  if (ier2 == 1 ) {
    stop(" Matrix with zero row in decompose" )
  } else if (ier2 == -1) {
    stop(" PARCOR for AR coefficients > 0.9" )
  } else if (ier2 == 2) {
    stop(" Singular matrix in decompose.zero divide in solve" )
  } else if (ier2 == 3) {
    stop(" Convergence in impruv.matrix is nearly singular" )
  } else if (ier2 ==400) {
    stop(" This model is unsuitable. The condition of the filtering routine is violated" )
  }

  trend <- NULL
  seasonal <- NULL
  ar <- NULL
  if (m3 == 0) arcoef <- NULL
  deff <- NULL

  sig2 <- 
  m12 <- m1 + m2*(period-1)
  xss <- array(z$xss, dim = c(mj, nmax))
  if (m1 > 0)
    trend <- xss[1,1:nmax]
  if (m2 > 0)
    seasonal <- xss[m1+1,1:nmax]
  if (m3 > 0)
    ar <- xss[m12+1, 1:nmax]
  if (m4 > 0)
    deff <- z$deff
  vss <- array(z$vss, dim = c(mj,mj,nmax))
  vss <- vss[1,1,1:nmax]

  yy <- y
  if (logt == 1) {
    yy <- log10(y)
  } else if (logt == 2) {
    yy <- log(y)
  }
  
  noise <- yy
  if (m1 > 0)
    noise <- noise - trend[1:n]
  if (m2 > 0)
    noise <- noise - seasonal[1:n]
  if (m3 > 0)
    noise <- noise - ar[1:n]
  if (m4 > 0)
    noise <- noise - deff[1:n]
  noise.list <- list(val = noise, range = c(ns, nfe), period = period)

  season.out <- list(tau2 = z$tau2[1:np], sigma2 = z$over, llkhood = z$ff,
                     aic = z$aic, trend = trend, seasonal = seasonal,
                     arcoef = arcoef, ar = ar, day.effect = deff,
                     noise = noise.list, cov = vss)
  class(season.out) <- "season"

  if (plot == TRUE) {
    rdata <- deparse(substitute(y))
    eval(parse(text=paste(rdata, "<- y")))
    if (logt == 1) {
      eval(parse(text=paste("plot.season(season.out, log10(", rdata, "), ...)")))
    } else if (logt == 2) {
      eval(parse(text=paste("plot.season(season.out, log(", rdata, "), ...)")))
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
  period <- x$noise$period

  if (is.null(rdata) == TRUE) {
    y <- NULL
    ts.atr <- c(1, 1, period)
    xtime <- "time"
  } else {
    y <- rdata
    tsname <- deparse(substitute(rdata))
    if (is.null(tsp(rdata)) == TRUE){
      y <- ts(y, start = 0, frequency = period)
    }
    ts.atr <- tsp(y)
    if (period == 4 || period == 12) {
      xtime <- "year"
    } else if (period == 24) {
      xtime <- "day"
    } else if (period == 5 || period == 7) {
      xtime <- "week"
    } else {
      xtime <- "time"
    }
  }

  trend <- x$trend
  seasonal <- x$seasonal
  ar <- x$ar
  day.effect <- x$day.effect
  noise <- x$noise$val
  noise <- ts(noise, start = ts.atr[1], frequency = ts.atr[3])

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

  old.par <- par(no.readonly = TRUE)
  par(mfrow = c(3, 1), xaxs = "i", yaxs = "i")
  nplot <- 0
  new.mar <- old.par$mar
  new.mar[1] <- new.mar[1] * 0.75
  new.mar[3] <- new.mar[3] * 0.5
  new.mgp <- old.par$mgp
  new.mgp[1] <- new.mgp[1] * 0.75

  nmax <- max(n, npe)
  xlim <- c(1, nmax)

### original and trend
  par(mar = new.mar, mgp = new.mgp)
  k <- NULL
  if (is.null(trend) == TRUE) {
    if (is.null(y) == FALSE) {
      ylim1 <- get.ylim(y)$minmax
      k <- get.ylim(y)$k
      mtitle <- paste(tsname)
      y <- ts(y, start = ts.atr[1], frequency = ts.atr[3])
      plot(y, type = "l", ylim = ylim1, xlab = xtime, ylab = "", main = mtitle, ...)
      nplot <- nplot + 1
    }

  } else {  # is.null(trend) == FALSE
    tt <- rep(NA, nmax)
    tt[ns:npe] <- trend[ns:npe]
    tt <- ts(tt, start = ts.atr[1], frequency = ts.atr[3])
    s2 <- x$sigma2
    t1 <- tt[1:nmax]
    t2 <- tt[1:nmax]
    for (i in 1:nmax) {
      t1[i] <- t1[i] - sqrt(x$cov[i] * s2)
      t2[i] <- t2[i] + sqrt(x$cov[i] * s2)
    }
    yrange <- range(t1[ns:nmax], t2[ns:nmax], y, na.rm = TRUE)
    ylim1 <- get.ylim(yrange)$minmax
    k <- get.ylim(yrange)$k

    trend <- ts(trend, start = ts.atr[1], frequency = ts.atr[3])
    if (is.null(y) == TRUE) {
      mtitle <- paste("Trend component")
      plot(trend, type = "n", ylim = ylim1, xlab = "", ylab = "",
           main = mtitle, ...)
      par(new = TRUE)
    } else {
      mtitle <- paste(tsname, "and trend component")
      y <- ts(y, start = ts.atr[1], frequency = ts.atr[3])
      yyy <- rep(NA, nmax)
      yyy[1:n] <- y[1:n]
      yyy <- ts(yyy, start = ts.atr[1], frequency = ts.atr[3])
      plot(yyy, type = "l", ylim = ylim1, xlab = "", ylab = "", main = mtitle, ...)
      par(new = TRUE)
    }
    plot(tt, type = "l", col = 2, 
         ylim = ylim1, xlab = xtime, ylab = "", main = "", xaxt = "n", ...) 
    par(new = TRUE)
    plot(t1, type = "l", lty = 3, col = 4,
         ylim = ylim1, xlab = "", ylab = "", main = "", xaxt = "n", ...)
    par(new = TRUE)
    plot(t2, type = "l", lty = 3, col = 4,
         ylim = ylim1, xlab = "", ylab = "", main = "", xaxt = "n", ...)
    nplot <- nplot + 1
  }
  par(mar = old.par$mar, mgp = old.par$mgp)


### seasonal
  if (is.null(seasonal) == FALSE) {
 
    if (is.null(k) == TRUE)
      k <- get.ylim(seasonal)$k
	  
    if (k == 0) {
	  ylim2 <- get.ylim(seasonal)$minmax * 2
	} else {
      ylim2 <- get.ylim(seasonal)$minmax
      if (ylim2[2] < 10 ** k) {
        ylim2[1] <- -10 ** k
        ylim2[2] <- 10 ** k
      }
    }

    seasonal <- ts(seasonal, start = ts.atr[1], frequency = ts.atr[3])
    ss <- rep(NA, nmax)
    ss[ns:npe] <- seasonal[ns:npe]
    ss <- ts(ss, start = ts.atr[1], frequency = ts.atr[3])
    plot(ss, type = "h", ylim = ylim2, xlab = xtime, ylab = "",
         main = "Seasonal component", ...) 
    nplot <- nplot + 1
  }


### AR
  if (is.null(ar) == FALSE) {

    if (is.null(k) == TRUE)
      k <- get.ylim(ar)$k
	  
    if (k == 0) {
      ylim3 <- get.ylim(ar)$minmax * 2
    } else {
      ylim3 <- get.ylim(ar)$minmax
      if (ylim3[2] < 10 ** k) {
        ylim3[1] <- -10 ** k
        ylim3[2] <- 10 ** k
      }
    }

    if (ylim3[1] > 0)
      ylim3[1] <- 0

    ar <- ts(ar, start = ts.atr[1], frequency = ts.atr[3])
    aa <- rep(NA, nmax)
    aa[ns:npe] <- ar[ns:npe]
    aa <- ts(aa, start = ts.atr[1], frequency = ts.atr[3])
    plot(aa, type = "h", xlab = xtime, ylim = ylim3, ylab = "",
         main = "AR component", ...)
    nplot <- nplot + 1
  }


### day.effect
  if (is.null(day.effect) == FALSE) {
  
    if (is.null(k) == TRUE)
      k <- get.ylim(day.effect)$k

    if (k == 0) {
      ylim4 <- get.ylim(day.effect)$minmax * 2
    } else {
      ylim4 <- get.ylim(day.effect)$minmax
      if (ylim4[2] < 10 ** k) {
        ylim4[1] <- -10 ** k
        ylim4[2] <- 10 ** k
      }
    }

    if (nplot == 3) {
      dev.new()
      par(mfrow = c(3, 1), xaxs = "i", yaxs = "i")
      nplot <- 0
    }
    day.effect <- ts(day.effect, start = ts.atr[1], frequency = ts.atr[3])
    dd <- rep(NA, nmax)
    dd[ns:npe] <- day.effect[ns:npe]
    dd <- ts(dd, start = ts.atr[1], frequency = ts.atr[3])
    plot(dd, type = "h", xlab = xtime, ylim = ylim4, ylab = "",
         main = "Trading day effect", ...)
    nplot <- nplot + 1
  }

### noise

  if (is.null(k) == TRUE)
    k <- get.ylim(noise)$k

  if (k == 0) {
    ylim5 <- get.ylim(noise)$minmax * 2
  } else {
    ylim5 <- get.ylim(noise)$minmax
    if (ylim5[2] < 10 ** k) {
      ylim5[1] <- -10 ** k
      ylim5[2] <- 10 ** k
    }
  }

  if (ylim5[1] > 0)
    ylim5[1] <- 0

  if (nplot == 3) {
    dev.new()
    par(mfrow = c(3, 1), xaxs = "i", yaxs = "i")
    nplot <- 0
  }

  sse <- rep(NA, nmax)
  sse[ns:nfe] <- noise[ns:nfe]
  sse <- ts(sse, start = ts.atr[1], frequency = ts.atr[3])

  plot(sse, type = "n", ylim = ylim5, xlab = xtime, ylab = "", main = "Noise", ...)
  par(new = TRUE)
  plot(sse, type = "h", xlab = "", ylim = ylim5, ylab = "", xaxt = "n", ...) 
  nplot <- nplot + 1

### Predicted values
  if (nfe < npe) {
    nps <- nfe + 1
    yy <- rep(NA, npe)
    for (i in nps:npe) {
      yy[i] <- 0
      if (is.null(trend) == FALSE)
        yy[i] <- yy[i] + trend[i]
      if (is.null(seasonal) == FALSE)
        yy[i] <- yy[i] + seasonal[i]
      if (is.null(ar) == FALSE)
        yy[i] <- yy[i] + ar[i]
      if (is.null(day.effect) == FALSE)
        yy[i] <- yy[i] + day.effect[i]
    }
    yy <- ts(yy, start = ts.atr[1], frequency = ts.atr[3])
    ydummy <- ts(rep(NA, nps), start=ts.atr[1], frequency = ts.atr[3])
    xlim <- c(0, nmax)

    if (nplot == 3) {
      dev.new()
      par(mfrow = c(3, 1), xaxs = "i", yaxs = "i")
    }
    par(mar = new.mar, mgp = new.mgp)

    if (is.null(trend) == FALSE) {
      t1 <- yy
      t2 <- yy
      for (i in nps:npe) {
        t1[i] <- t1[i] - sqrt(x$cov[i] * s2)
        t2[i] <- t2[i] + sqrt(x$cov[i] * s2)
      }
      ylim6 <- range(ylim1, t1[nps:npe], t2[nps:npe], na.rm = TRUE)

      plot(yy, type = "l", col = 2, ylim = ylim6,
           xlab = xtime, ylab = "", main = "Predicted values", xaxt = "n", ...)
      par(new = TRUE)

      if (is.null(y) == FALSE) {
        y0 <- rep(NA, npe)
#        y0[ns:nfe] <- y[ns:nfe]
        y0[1:length(y)] <- y
        y0 <- ts(y0, start = ts.atr[1], frequency = ts.atr[3])
        plot(y0, type = "l", ylim = ylim6, xlab = "", ylab = "", main = "", ...)
        par(new = TRUE)
      }

      plot(t1, type = "l", col = 4, ylim = ylim6,
           xlab = "", ylab = "", main = "", xaxt = "n", ...)
      par(new = TRUE)
      plot(t2, type = "l", col = 4, ylim = ylim6,
           xlab = "", ylab = "", main = "", xaxt = "n", ...)
      abline(v = tsp(ydummy)[2], lty = 2)
	  
    } else {
      ylim6 <- range(ylim1, yy[nps:npe], na.rm = TRUE)
      if (is.null(y) == FALSE) {
        y0 <- rep(NA, npe)
        y0[ns:nfe] <- y[ns:nfe]
        y0 <- ts(y0, start = ts.atr[1], frequency = ts.atr[3])
        plot(y0, type = "l", ylim = ylim6, xlab = "", ylab = "", main = "", ...) 
        par(new = TRUE)
      }
      plot(yy, type = "l", col = 2, ylim = ylim6, xaxt = "n",
           xlab = xtime, ylab = "", main = "Predicted values", ...)
      abline(v = tsp(ydummy)[2], lty = 2)
    }
    par(mar = old.par$mar, mgp = old.par$mgp)
  }
  par(mfrow = old.par$mfrow, xaxs = old.par$xaxs)
}


get.ylim <- function(y) {
  ymax <- max(y)
  ymin <- min(y)
  yd <- ymax - ymin
  k <- floor(log10(yd))
  ymax <- ((ymax %/% 10^k) + 1) * (10 ** k)
  ymin <- (ymin %/% 10^k) * (10 ** k)
  if (ymin < 0) {
    ys <- max(abs(ymin), abs(ymax))
    ymin <- -ys
    ymax <- ys
  }
  return(list(minmax = c(ymin, ymax), k = k))
}
