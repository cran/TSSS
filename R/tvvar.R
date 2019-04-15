# PROGRAM 13.1  TVVAR
tvvar <- function(y, trend.order, tau2.ini = NULL, delta, plot = TRUE, ...)
{
  n <- length(y)             # length of data
  m <- trend.order           # trend order
# tau2.ini                   # initial variance of systen noise
# delta                      # delta for computing variance of system noise
  iopt <- 1                  # search method
  if (is.null(tau2.ini)) {
    iopt  <- 0
    tau2.ini <- 0
    delta <- 0
  }

  z <- .Call("TvvarC",
             as.double(y),
             as.integer(n),
             as.integer(m),
             as.double(tau2.ini),
             as.integer(iopt),
             as.double(delta))

  nordata <- z[[2L]]
  sm <- z[[3L]]
  trend <- array(z[[5L]], dim = c(z[[4L]], 3))
  noise <- z[[6L]]

  tvvar.out <- list(tvv = z[[1L]], nordata = nordata, sm = sm, trend = trend,
                    noise = noise, tau2 = z[[7L]], sigma2 = z[[8L]],
                    llkhood = z[[9L]], aic = z[[10L]],
                    tsname = deparse(substitute(y)))
  class(tvvar.out) <- "tvvar"

  if (plot)
    plot.tvvar(tvvar.out, ...)
  return(tvvar.out)
}

print.tvvar <- function(x, ...)
{
  message(gettextf("\n tau2\t\t%12.5e", x$tau2), domain = NA)
  message(gettextf(" sigma2\t\t%12.5e", x$sigma2), domain = NA)
  message(gettextf(" log-likelihood\t%12.3f", x$llkhood), domain = NA)
  message(gettextf(" aic\t\t%12.3f\n", x$aic), domain = NA)
}

plot.tvvar <- function(x, ...)
{
    old.par <- par(no.readonly = TRUE)
##    par(mfrow=c(4,1))
    par(mfrow=c(3,1), xaxs = "i")

##    plot(nordata, type = "h", xlab = "m", ylab = "normalized data")
##    plot(nordata, type = "l", xlab = "m", ylab = "normalized data")
##    par(new = TRUE)
##    abline(h = 0)

##    plot(sm, type = "l", xlab = "m", ylab = "s(m)=y(2m-1)**2 + y(2m)**2",
    plot(x$sm, type = "l", xlab = "m", main = paste(x$tsname),
         ylab = "transformed data", ...)

    ylim <- c(floor(min(x$trend)), ceiling(max(x$trend)))
    plot(x$trend[, 1], type = "l", ylim = ylim, xlab = "m", ylab = "trend", ...) 
    par(new = TRUE) 
    plot(x$trend[, 2], type = "l", col=2, ylim = ylim, xlab = "", ylab = "", ...)
    par(new = TRUE)
    plot(x$trend[, 3], type = "l", ylim = ylim, xlab = "", ylab = "", ...) 

##    plot(noise, type = "h", xlab = "m", ylab = "noise")
    plot(x$noise, type = "h", xlab = "m", ylab = "standardized data ", ...)
    par(new = TRUE)
    abline(h = 0)
    par(old.par)
}

# PROGRAM 13.2  TVAR
tvar <- function (y, trend.order = 2, ar.order = 2, span, outlier = NULL,
                  tau2.ini = NULL, delta, plot = TRUE)     
{
# y                     # original data
  n <- length(y)        # data length
# ar.order              # AR order
# trend.order           # Trend order
  if (trend.order!=1 && trend.order!=2)
    stop("'trend.order' is 1 or 2" )
# span                  # local stationary span
# outlier               # position of i-th outlier
  if (is.null(outlier)) {
    nout <- 0           # number of outliers
    outlier <- 0
  } else {
    nout <- length(outlier)
  }
# tau2.ini              # initial variance of systen noise
# delta                 # delta for computing variance of system noise
  iopt <- 1             # search method
  if (is.null(tau2.ini)) {
    iopt  <- 0
    tau2.ini <- 0
    delta <- 0
  }

  z <- .Call("TvarC",
             as.double(y),
             as.integer(n),
             as.integer(ar.order),
             as.integer(trend.order),
             as.integer(span),
             as.integer(iopt),
             as.integer(nout),
             as.integer(outlier),
             as.double(tau2.ini),
             as.double(delta))

  nn <- n/span
  sigma2 <- z[[2L]]
  arcoef <- array(z[[5L]], dim = c(ar.order, nn))
  parcor <- array(z[[6L]], dim = c(ar.order, nn))
  xx <- span
  for (i in 2:nn)
    xx <- c(xx, i*span)

  tvar.out <- list(arcoef = arcoef, sigma2 = sigma2, tau2 = z[[1L]],
                   llkhood = z[[3L]], aic = z[[4L]], parcor = parcor)
  class(tvar.out) <- "tvar"

  if (plot)
    plot.parcor(parcor, span)

  return(tvar.out)
}

print.tvar <- function(x, ...)
{
  message(gettextf("\n tau2\t\t%12.5e", x$tau2), domain = NA)
  message(gettextf(" sigma2\t\t%12.5e", x$sigma2), domain = NA)
  message(gettextf(" log-likelihood\t%12.3f", x$llkhood), domain = NA)
  message(gettextf(" aic\t\t%12.3f\n", x$aic), domain = NA)
}

plot.parcor <- function(parcor, span)
{
  ar.order <- dim(parcor)[1]
  n <- dim(parcor)[2]

  x <- span
  for (i in 2:n)
    x <- c(x, i*span)

  if (ar.order > 3)
    nc <- 3
  old.par <- par(no.readonly = TRUE)
  par(mfrow = c(nc, 1))

  for (i in 1:ar.order ) {
    if ((i%%3 == 1) & i > 1)
      par(ask = TRUE)
      plot(x, parcor[i, ], type = "l", xlab = "", ylim = c(-1.0, 1.0),
           ylab = paste("parcor[", i, ", ]"), xaxs = "i", yaxs = "i")
  }

  par(old.par)
}


# PROGRAM 13.3  TVSPC
tvspc <- function (arcoef, sigma2, var = NULL, span = 20, nf = 200)
{
# arcoef                        # Time varying AR coefficient
  n <- dim(arcoef)[2]           # number of points
  ar.order <- dim(arcoef)[1]    # AR order

  if (is.null(var)) {
    ivar <- 0
    var <- rep(0, n * span)     # time varying variance
  } else {
    ivar <- 1                   # =1: for variance correction
  }
# span                          # local stationary span
# nf                            # number of frequencies

  z1 <- .Call("TvspcC",
              as.integer(n),
              as.integer(ar.order),
              as.integer(span),
              as.integer(nf),
              as.integer(ivar),
              as.double(sigma2),
              as.double(arcoef),
              as.double(var))

  spectra <- array(z1[[1L]], dim = c(nf+1, n))
  yy <- seq(0, n*span, length = n)

  out <- list(y = yy, z = spectra)
  class(out) <- c("tvspc")
  return(out)

}

# PROGRAM 13.3  TVSPC
plot.tvspc <- function (x, theta = 0, phi = 15, expand = 1, col = "lightblue",
                        ticktype= "detail", ...)
{
#  spectra <- x$z
#  nf1 <- dim(spectra)[1]
  nf1 <- dim(x$z)[1]
  xx <- seq(0, 0.5, length = nf1)
#  zlim <- c(floor(spectra), ceiling(spectra))

  persp(xx, x$y, x$z, xlab = "f", ylab = "n", zlab = "log p(f)", theta = theta,
        phi = phi, expand = expand, col = col, ticktype = ticktype, ...)
}

