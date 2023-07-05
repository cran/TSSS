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

  n2 <- n / 2
  
  z <- .Fortran(C_tvvarf,
             as.double(y),
             as.integer(n),
             as.integer(m),
             as.double(tau2.ini),
             as.integer(iopt),
             as.double(delta),
             tvvar = double(n2),
             nordat = double(n),
             y1 = double(n2),
             n1 = integer(1),
             trend = double(n2 * 3),
             noise = double(n2),
             taumax = double(1),
             sig2m = double(1),
             ffmax = double(1),
             aic = double(1))

  trend <- array(z$trend, dim = c(z$n1, 3))

  tvvar.out <- list(tvv = z$tvvar, nordata = z$nordat, sm = z$y1, trend = trend,
                    noise = z$noise, tau2 = z$taumax, sigma2 = z$sig2m,
                    llkhood = z$ffmax, aic = z$aic,
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

    new.mar <- old.par$mar
    new.mar[1] <- 3.1
    new.mar[3] <- 1.1
    new.mgp <- old.par$mgp
    new.mgp[1] <- 2.25
    new.mgp[2] <- 0.75

    par(oma = c(0, 0, 2, 0), mfrow=c(4, 1), xaxs = "i", mar = new.mar,
        mgp = new.mgp)

    ylim <- c(floor(min(x$trend, x$sm)), ceiling(max(x$trend, x$sm)))
    plot(x$sm, type = "l", ylim = ylim, xlab = "", main = "",
         ylab = "transformed data", mar = new.mar)
    title(xlab = "m", mgp = c(2, 0.75, 0))

    plot(x$trend[, 1], type = "l", ylim = ylim, xlab = "", ylab = "trend t(m)") 
    par(new = TRUE) 
    plot(x$trend[, 2], type = "l", col=2, ylim = ylim, xlab = "", ylab = "")
    par(new = TRUE)
    plot(x$trend[, 3], type = "l", ylim = ylim, xlab = "", ylab = "") 
    title(xlab = "m", mgp = c(2, 0.75, 0))

    plot(x$noise, type = "h", xlab = "", ylab = "residual")
    par(new = TRUE)
    abline(h = 0)
    title(xlab = "m", mgp = c(2, 0.75, 0))

    ylim <- c(floor(min(x$nordat)), ceiling(max(x$nordat)))
    plot(x$nordat, type = 'h', ylim = ylim, xlab = "", ylab = "normalized data")
    abline(h=0)
    title(xlab = "m", mgp = c(2, 0.75, 0))

    mtext(paste(x$tsname), side = 3, line = 0, outer = TRUE)

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
  if (trend.order != 1 && trend.order != 2)
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

  nn <- n/span
  
  z <- .Fortran(C_tvar,
                as.double(y),
                as.integer(n),
                as.integer(ar.order),
                as.integer(trend.order),
                as.integer(span),
                as.integer(iopt),
                as.integer(nout),
                as.integer(outlier),
                as.double(tau2.ini),
                as.double(delta),
                tau2 = double(1),
                sigma2 = double(1),
                lkhood = double(1),
                aic = double(1),
                arcoef = double(nn * ar.order),
                parcor = double(nn * ar.order))

  arcoef <- array(z$arcoef, dim = c(ar.order, nn))
  parcor <- array(z$parcor, dim = c(ar.order, nn))

  tvar.out <- list(arcoef = arcoef, sigma2 = z$sigma2, tau2 = z$tau2,
                   llkhood = z$lkhood, aic = z$aic, parcor = parcor)
  class(tvar.out) <- "tvvar"

  if (plot)
    plotParcor(parcor, span)

  return(tvar.out)
}

plotParcor <- function(parcor, span)
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

  nf1 <- nf + 1
  
  z <- .Fortran(C_tvspc,
                 as.integer(n),
                 as.integer(ar.order),
                 as.integer(span),
                 as.integer(nf),
                 as.integer(ivar),
                 as.double(sigma2),
                 as.double(arcoef),
                 as.double(var),
                 spec = double(nf1 * n))

  spectra <- array(z$spec, dim = c(nf1, n))
  yy <- seq(0, n*span, length = n)

  out <- list(y = yy, z = spectra)
  class(out) <- c("tvspc")
  return(out)

}

# PROGRAM 13.3  TVSPC
plot.tvspc <- function(x, tvv = NULL, dx = 2, dy = 0.25, ...) {

  zoff <- FALSE

  if (dx < 0 || dx > 10)
    stop("0 <= dx <= 10")
  if (dy > 1 || dy < 0.05)
    stop("0.05 <= dy <= 1")

  nf1 <- dim(x$z)[1]
  nf <- nf1 - 1
  if (zoff == FALSE) {
    xoff <- 0
  } else {
    xoff <- nf/20
  }

  nrep <- dim(x$z)[2]
  zmax <- max(x$y)

  span <- zmax / nrep
  xw <- (nf + dx * nrep)
  xmax <- max(xw, (2.5*nf)) * 1.1

  z <- x$z
  if (is.null(tvv) == FALSE) {
    if ((length(tvv) /10) != nrep)
      stop("the length of time varying variance is invalid.")
    for (i in 1:nrep)
      z[, i] <- z[, i] + log10(tvv[i * 10])
  }
  yoff <- 1
  ymin <- trunc(min(z)) - yoff
  ymax0 <- ceiling(max(z))

  yw <- ymax0 + dy * nrep - ymin
  ymax <- yw * 1.1

  iexp <- 0
  zmax0 <- zmax
  while (zmax0 > 10) {
    zmax0 <- zmax0 / 10
    iexp <- iexp + 1
  }
  dz <- 10 ** iexp
  if (zmax0 > 8) {
    dz <- dz * 2
  } else if (zmax0 < 3) {
    while (zmax0 < 3) {
      dz <- dz / 2
      zmax0 <- zmax0 * 2
    }
  }

  x <- c(0:nf)
  y <- z[, 1]
  zz <- y

  old.par <- par(no.readonly = TRUE)
  par(bty = "n", xaxs = "i", yaxs = "i")
  new.mgp <- c(1.8, 0.5, 0)

# plot the spectrum
  plot(x, zz, type="l", xlim = c(0, xmax), ylim = c(ymin, ymax), xaxt = "n",
       yaxt = "n", xlab = "", ylab = "", ...)
  title(xlab = "f", adj = nf / 2 / xmax, mgp = new.mgp)
  title(ylab = "log p(f)", adj = (ymax0 - ymin) / 2 / ymax, , mgp = new.mgp)

  for (i in 2:nrep){
    par(new=T)
    x <- x + dx
    for (j in 1:nf1)
      zz[j] <- z[j,i] + dy*(i-1)
    for (j in 1:(nf-dx+1))
      zz[j] <- max(zz[j], y[j+dx])
    y <- zz
    plot(x, zz, type="l", xlim = c(0, xmax), ylim = c(ymin, ymax), xaxt = 'n',
         yaxt = "n", xlab = "", ylab = "", ...)
  }
  xend <- x[nf1]
  yend <- dy * nrep

# draw axes
  xlab <- c("0", "","","","", "0.5")
  axis(1, at = seq(0, nf, length = 6), labels = paste(xlab), tck = -0.01,
       cex.axis = 0.8, mgp = new.mgp)

  ylab <- c(ymin:ymax0)
  axis(2, at = ymin:ymax0, cex.axis = 0.8, mgp = new.mgp)

  zlab <- seq(0, zmax, by = dz)
  nz <- length(zlab)
  zlab.max <- zlab[nz]
  maxp <- zlab.max / span
  zlab <- as.character(zlab)
  segments(x0 = nf + xoff, y0 = ymin, x1 = xend + xoff, y1 = ymin + yend)

  ddx <- (dx * maxp) / (nz - 1)
  ddy <- (dy * maxp)  / (nz - 1)

  for (i in 1:nz) {
    xx <- nf + ddx * (i-1) + xoff
    yy <- ymin + ddy *(i-1)
    if (i != 1)
      text(xx + 5, yy, zlab[i], cex=0.8, pos=4)
    segments(xx, yy, xx + 5, yy, xpd = TRUE) 
  }
  text(xend + xoff + 5 , max(yend, yy+ddy/3), "n", cex = 0.8, pos = 4)

  par(bty = old.par$bty, xaxs = old.par$xaxs, yaxs = old.par$yaxs)
}
