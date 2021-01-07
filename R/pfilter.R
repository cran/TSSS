pfilter <- function(y, m = 10000, model = 0, lag = 20, initd = 0, sigma2, tau2,
                    alpha = 0.99, bigtau2 = NULL, init.sigma2 = 1,
                    xrange = NULL, seed = NULL, plot = TRUE,  ...)
{
# y                     # original data
  n <- length(y)        # data length
# m                     # number of particles
# model                 # model (0,1,2)
  if (model < 0 || model > 2)
    stop("'model' must be 0, 1 or 2.")
# lag                   # fixed lag smoother
# initd                 # initial distribution (0,1,2,3)
  if (initd < 0 || initd > 3)
    stop("'model' must be 0, 1, 2 or 3.")
# sigma2                # observation noise variance
  if (sigma2 < 0 || sigma2 == 0)
    stop("'sigma2' must be greater than zero.")
# tau2                  # system noise variance
  if (tau2 < 0)
    stop("'tau2' must be greater than or equal to zero.")
# alpha                 # mixture weight (model = 2)
# xrange                # lower and upper bounds of the distribution's range
  if (is.null(xrange) == TRUE) {
    r2 <- floor((max(y) - min(y)) /  2)
    xmin <- floor(min(y) - r2)
    xmax <- ceiling(max(y) + r2)
  } else {
    xmin <- xrange[1]
    xmax <- xrange[2]
    if (xmin == xmax || xmin > xmax)
      stop("xrange[2] must be greater than xrange[1].")
  }
# bigtau2               # vairance of the second component (model = 2)
  if (is.null(bigtau2) == TRUE) {
    bigtau2 <- (xmax - xmin) / 2
  } else if (bigtau2 < 0) {
    stop("'bigtau2' must be greater than or equal to zero.")
  }
# init.sigma2           # variance of initial state distribution (initd = 0)
  if (init.sigma2 < 0)
    stop("'init.sigma2' must be greater than or equal to zero.")
# seed :     for random number generator
  ix <- seed
  if (is.null(ix))
    ix <- -1

  z <- .Call("PfiltC",
             as.double(y),
             as.integer(n),
             as.integer(m),
             as.integer(model),
             as.integer(lag),
             as.integer(initd),
             as.double(sigma2),
             as.double(tau2),
             as.double(alpha),
             as.double(bigtau2),
             as.double(init.sigma2),
             as.double(xmin),
             as.double(xmax),
             as.integer(ix))

  t <- array(z[[1L]], dim = c(n, 8))
  tdist <- t[, 7:1]
  ff <- z[[2L]]
  message(gettextf("\n log-likelihood\t%12.3f", ff), domain = NA)

  pfilter.out <- list(llkhood = ff, smooth.dist = tdist)
  class(pfilter.out) <- "pfilter"

  if (plot) {
    plot.pfilter(pfilter.out, ...)
    invisible(pfilter.out)
  } else {
    return(pfilter.out)
  }
}


pfilterNL <- function(y, m = 10000, lag = 20, sigma2, tau2, xrange = NULL,
                      seed = NULL, plot = TRUE,  ...)
{
# y                     # original data
  n <- length(y)        # data length
# m                     # number of particles
# lag                   # fixed lag smoother
# sigma2                # observation noise variance
  if (sigma2 < 0 || sigma2 == 0)
    stop("'sigma2' must be greater than zero.")
# tau2                  # system noise variance
  if (tau2 < 0)
    stop("'tau2' must be greater than or equal to zero.")
# xrange                # lower and upper bounds of the distribution's range
  if (is.null(xrange) == TRUE) {
    r2 <- floor((max(y) - min(y)) /  2)
    xmin <- floor(min(y) - r2)
    xmax <- ceiling(max(y) + r2)
  } else {
    xmin <- xrange[1]
    xmax <- xrange[2]
    if (xmin == xmax || xmin > xmax)
      stop("xrange[2] must be greater than xrange[1].")
  }
# seed :     for random number generator
  ix <- seed
  if (is.null(ix))
    ix <- -1

  z <- .Call("PfiltnlC",
             as.double(y),
             as.integer(n),
             as.integer(m),
             as.integer(lag),
             as.double(sigma2),
             as.double(tau2),
             as.double(xmin),
             as.double(xmax),
             as.integer(ix))

  t <- array(z[[1L]], dim = c(n, 8))
  tdist <- t[, 7:1]
  ff <- z[[2L]]
  message(gettextf("\n log-likelihood\t%12.3f", ff), domain = NA)

  pfilterNL.out <- list(llkhood = ff, smooth.dist = tdist)
  class(pfilterNL.out) <- "pfilter"

  if (plot) {
    plot.pfilter(pfilterNL.out, ...)
    invisible(pfilterNL.out)
  } else {
    return(pfilterNL.out)
  }
}


plot.pfilter <- function(x, ...)
{
  t <- x$smooth.dist
  mint <- floor(min(t) - 0.2)
  maxt <- ceiling(max(t) + 0.2)
#  minmaxt <- max(abs(mint), abs(maxt))
#  ylim <- c(-minmaxt, minmaxt)
  ylim <- c(mint, maxt)
  n <- dim(t)[1]
  xx <- seq(1,n, length=n)
  gcol <- c("gray90", "gray80", "gray70", "gray70", "gray80", "gray90")
  mtitle <- "Marginal smoothed distribution of trend"

  plot(t[, 1], type = 'l', ylim = ylim, main = mtitle, xlab = "time",
       ylab = "x", xaxs = 'i', yaxs = 'i', ...)
  for(i in 2:7) {
    par(new = TRUE)
    plot(t[, i], type = 'l', ylim = ylim, xlab = "", ylab = "", xaxs = 'i',
         yaxs = 'i', ...)
    y1 <- t[xx,i-1]
    y2 <- t[xx,i]
    polygon( c(xx,rev(xx)), c(y1,rev(y2)), col=gcol[i-1])
  }
  par(new = TRUE)
  plot(t[, 4], type = 'l', ylim = ylim, xlab = "", ylab = "", xaxs = 'i',
       yaxs = 'i', lwd = 2, ...)
}

