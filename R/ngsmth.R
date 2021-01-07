# PROGRAM 14.1  NGSMTH
ngsmth <- function(y, noisev = 2,tau2, bv = 1.0, noisew = 1, sigma2, bw = 1.0, 
                   initd = 1, k = 200, plot = TRUE, ...)
{
# y                     # original data
  n <- length(y)        # data length
# noisev                # type of system noise density (0,1,2,3)
                        # 1: Gaussian (normal) / 2: Pearson family
                        # 3: two-sides exponential
  if (noisev!=1 && noisev!=2 && noisev!=3)
    stop("'noisev' is numeric in {1,2,3}" )
# tau2                  # variance of dispersion of system noise
# bv                    # shape parameter of system noise (for noisev=2)
  if (noisev == 2)
    if ((bv == 0.5) || (bv < 0.5))
      stop("shape parameter of system noise is greater than 0.5")
# noisew                # type of observation noise density (0,1,2,3,4)
                        # 1: Gaussian (normal) / 2: Pearson family
                        # 3: two-sided exponential / 4:double exponential
  if (noisew!=1 && noisew!=2 && noisew!=3 && noisew!=4)
    stop("'noisew' is numeric in {1,2,3,4}" )
# sigma2                # variance of dispersion of observation noise
# bw                    # shape parameter of observation noise (for noisew=2)
  if (noisew == 2)
    if ((bw == 0.5) || (bw < 0.5))
      stop("shape parameter of observation noise is greater than 0.5")
# initd                 # type of density function
                        # 1: Gaussian (normal) / 2: uniform
                        # 3: two-sided exponential
  if (initd!=1 && initd!=2 && initd!=3)
    stop("'initd' is numeric in {1,2,3}")
  if (initd==3)
    initd <- 0
# k                     # number of intervals

  ns <- 1
  nfe <- n
  npe <- n
  k1 <- k+1

  z <- .Call("NgsmthC",
             as.double(y),
             as.integer(n),
             as.integer(noisev),
             as.double(tau2),
             as.double(bv),
             as.integer(noisew),
             as.double(sigma2),
             as.double(bw),
             as.integer(initd),
             as.integer(ns),
             as.integer(nfe),
             as.integer(npe),
             as.integer(k1))

  trend <- array(z[[1L]], dim = c(n, 7))
  smt <- array(z[[2L]], dim = c(k1, n))
  ll <- z[[3L]]
  message(gettextf("\n log-likelihood\t%12.3f", ll), domain = NA)

  ngsmth.out <- list(trend = trend, smt = smt)
  class(ngsmth.out) <- "ngsmth"

  if (plot)
    plot.ngsmth(ngsmth.out, "trend",...)

  return(ngsmth.out)
}

plot.ngsmth <- function(x, type = c("trend", "smt"), theta = 0, phi = 15,
                        expand = 1, col = "lightblue", ticktype= "detail", ...)
{
  ntype <- length(type)

  if (ntype == 0 )
    return()
  type1 <- FALSE
  type2 <- FALSE
  for (i in 1:ntype) {
    if (type[i] == "trend")
      type1 <- TRUE
    if (type[i] == "smt")
      type2 <- TRUE
  }

  if (type1 == TRUE) {
    old.par <- par(no.readonly = TRUE)
    ylim <- c(floor(min(x$trend)), ceiling(max(x$trend)))
    for (i in 1:7) {
      if (i != 1)
        par(new = TRUE) 
      if (i == 4) {
        plot(x$trend[, i], type = "l", col = 2, ylim = ylim, xlab = "n",
             ylab = "trend  tn", xaxs = "i", yaxs = "i", ...)
      } else {
        plot(x$trend[, i], type = "l", ylim = ylim, xlab = "", ylab = "",
             xaxs = "i", yaxs = "i", ...)
      }
    }
    par(old.par)
  }

  if (type2 == TRUE) {
    if (type1 == TRUE)
      dev.new()

    k <- dim(x$smt)[1] - 1
    n <- dim(x$smt)[2]

    ndif <- 1
    if (n >= 100)
      ndif <- 2
    if (n >= 200)
      ndif <- 4
    if (n >= 300)
      ndif <- 6
    if (n >= 500)
      ndif <- as.integer(n / 50)
    n0 <- as.integer(ndif / 2) + 1
    nn <- as.integer((n-n0) / ndif) + 1
    ss <- array(0, dim = c(k+1, nn))
    jj <- 0
    for (j in n0:n) {
      if ((j-n0)%%ndif == 0) {
        ss[, jj] <- x$smt[, j]
        jj <- jj + 1
      }
    }

    xs <- floor(min(x$trend))
    xe <- ceiling(max(x$trend))
    xx <- seq(xs, xe, length = k+1)
    yy <- seq(1, n, length = nn)
    zlim <- c(0, ceiling(max(ss)))
    persp(xx, yy, ss, zlim = zlim, xlab = "tn", ylab = "n", zlab = "p(tn)",
          theta = theta, phi = phi, expand = expand, col = col,
          ticktype = ticktype, ...)
  }  
}
