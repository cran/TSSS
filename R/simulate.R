# PROGRAM 15.1
simssm <- function(n = 200, trend = NULL, seasonal.order = 0, seasonal = NULL, arcoef = NULL,
                   ar = NULL, tau1 = NULL, tau2 = NULL, tau3 = NULL, 
                   sigma2 = 1.0, seed = NULL, plot = TRUE, ...)
{
  k <- 0
  if (is.null(trend)) {
    m1 <- 0
    x <- NULL
    tau1 <- 0
  } else {
    m1 <- length(trend)      # trend order   1 or 2
    if (is.null(tau1))
      stop("'tau1' is not specified")
    if (m1 > 2)
      stop("trend order is 0, 1 or 2")
    x <- trend
    k <- k + 1
  }

  m2 <- seasonal.order 
  if (m2 == 0) {
    period <- 0
    tau2 <- 0
  } else if (m2 == 1 || m2 == 2) {
    if (is.null(seasonal))
      stop("'seasonal' is not specified")
    if (is.null(tau2))
      stop("'tau2' is not specified")
    period <- length(seasonal) + 1
    x <- c(x, seasonal)
    k <- k + 1
  } else {
    stop("seasonal.order is 0, 1 or 2")
  }

  m3 <- length(arcoef)
  if (m3 == 0) {
    ar <- 0
    tau3 <- 0
  } else {
    if (is.null(tau3))
      stop("'tau3' is not specified")
    if (is.null(ar)) {
      ar <- rep(0, m3)
    } else {
      if(length(ar) != m3)
        stop("'ar' is invalid")
    }
    x <- c(x, ar)
    k <- k + 1
  }

  m <- m1 + m2 * (period - 1) + m3

# n:    simulation interval
# seed :     for random number generator
  ix <- seed
  if (is.null(ix))
    ix <- -1

  z <- .Call("SimssmC",
             as.integer(m1),
             as.integer(m2),
             as.integer(m3),
             as.integer(m),
             as.integer(k),
             as.integer(n),
             as.integer(ix),
             as.double(sigma2),
             as.integer(period),
             as.double(tau1),
             as.double(tau2),
             as.double(tau3),
             as.double(arcoef),
             as.double(x))

  out <- z[[1L]]
  class(out) <- "simulate"

  if (plot) {
    plot.simulate(out, c(1, n), ...)
    invisible(out)
  } else out
}


# PROGRAM 15.2
ngsim <- function(n = 200, trend = NULL, seasonal.order = 0, seasonal = NULL, arcoef = NULL,
                  ar = NULL, noisew = 1, wminmax = NULL, paramw = NULL,
                  noisev = 1, vminmax = NULL, paramv = NULL, seed = NULL,
                  plot = TRUE, ... )
{
  k <- 0
  if (is.null(trend)) {
    m1 <- 0
    x <- NULL
  } else {
    m1 <- length(trend)     # trend order   1 or 2
    if (m1 > 2)
      stop("trend order is 0, 1 or 2")
    x <- trend
    k <- k + 1
  }

  m2 <- seasonal.order      # seasonal order
  if (m2 == 0) {
    period <- 0
  } else if (m2 == 1 || m2 == 2) { 
    if (is.null(seasonal))
      stop("'seasonal' is not specified") 
    period <- length(seasonal) + 1
    x <- c(x, seasonal)
    k <- k + 1
  } else {
    stop("seasonal.order is 0, 1 or 2")
  }

  m3 <- length(arcoef)
  if (m3 == 0) {
    ar <- 0
  } else {
    if (is.null(ar)) {
      ar <- rep(0, m3)
    } else {
      if(length(ar) != m3)
        stop("'ar' is invalid")
    }
    x <- c(x, ar)
    k <- k + 1
  }

  m <- m1 + m2 * (period - 1) + m3

# n:    simulation interval
# seed :     for random number generator
  ix <- seed
  if (is.null(ix))
    ix <- -1

# noisew:    type of observational noise
  if (noisew < -3 || noisew > 3)
    stop("type of observational noise is one of -3, -2, -1, 0, 1, 2, 3")
# noisev:    type of system noise
  if (noisev < -3 || noisev > 3)
    stop("type of system noise is one of -3, -2, -1, 0, 1, 2, 3")

# wparam:      parameter of the observational noise density
# vparam:      parameter of the system noise density
  pw <- c(0.0, 1.0, 1.0)
  pv <- c(0.0, 1.0, 1.0)
    
  if (noisew == 1) {
    if (is.null(paramw))
      stop("'paramw' is not specified")
    if (paramw[1] < 0)
      stop("the value of 'paramw' (variance) is non-negative")
    pw[2] <- paramw[1]
  } else if (noisew == 2) {
    if (is.null(paramw))
      stop("'paramw' is not specified")
    if (length(paramw) != 2)
      stop("'paramw' is not specified")
    if ((paramw[1] == 0) || (paramw[1] < 0))
      stop("dispersion parameter of the observational noise density is greater
            than 0")
    if ((paramw[2] == 0) || (paramw[2] < 0))
      stop("shape parameter of the observational noise density is greater than
            0")
    pw[2:3] <- paramw[1:2]
  }

  if (noisev == 1) {
    if (is.null(paramv))
      stop("'paramv' is not specified")
    if (paramv[1] < 0)
      stop("the value of 'paramv' (variance) is non-negative")
    pv[2] <- paramv[1]
  } else if (noisev == 2) {
    if (is.null(paramv))
      stop("'paramv' is not specified")
    if (length(paramv) != 2)
      stop("'paramv' is not specified")
    if ((paramv[1] == 0) || (paramv[1] < 0))
      stop("dispersion parameter of the system noise density is greater than 0")
    if ((paramv[2] == 0.5) || (paramv[2] < 0.5))
      stop("shape parameter of the system noise density is greater than 0.5")
    pv[2:3] <- paramv[1:2]
  }

# wminmax:     lower and upper bound of observational  noise
  if (is.null(wminmax)) {
    if (pw[3] > 1.6) {
      sd4 <- 4 * sqrt(pw[2]**2 / (2*pw[3] - 3))
      wmin <- -sd4
      wmax <- sd4
    } else {
      wmin <- -20
      wmax <- 20
    }
  } else {
    wmin <- wminmax[1]
    wmax <- wminmax[2]
  }
# vminmax:     lower and upper bound of system noise
  if (is.null(vminmax)) {
    if (pv[3] > 1.6) {
      sd4 <- 4 * sqrt(pv[2]**2 / (2*pv[3] - 3))
      vmin <- -sd4
      vmax <- sd4
    } else {
      vmin <- -20
      vmax <- 20
    }
  } else {
    vmin <- vminmax[1]
    vmax <- vminmax[2]
  }

  z <- .Call("NgsimC",
             as.integer(m1),
             as.integer(m2),
             as.integer(m3),
             as.integer(m),
             as.integer(k),
             as.integer(n),
             as.integer(ix),
             as.integer(noisew),
             as.double(wmin),
             as.double(wmax),
             as.double(pw),
             as.integer(noisev),
             as.double(vmin),
             as.double(vmax),
             as.double(pv),
             as.integer(period),
             as.double(arcoef),
             as.double(x))

  out <- z[[1L]]
  class(out) <- "simulate"

  if (plot) {
    plot.simulate(out, c(1, n), ...)
    invisible(out)
  } else out

}

plot.simulate <- function(x, use = NULL, ...)
{
  n <- length(x)
  if (is.null(use) == TRUE) {
    st <- 1
    ed <- n
  } else {
    st <- use[1]
    ed <- use[2]
  }
  if (st < 1)
      stop("start time must be between 1 and n.\n")
  if (ed > n)
      stop("end time must be between 1 and n.\n")
  if (ed > n)
      stop("used time interval is invalid.\n")

  xx <- c(st:ed)
  yy <- x[st:ed]
  plot(xx, yy, type = "l", xlab = "", ylab = "y(n)",
       main = "Simulated time series", xaxs = "i", ...)
}
