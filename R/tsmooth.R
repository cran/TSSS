# PROGRAM 9.1  SMOOTH
tsmooth <- function( y, f, g, h, q, r, x0 = NULL, v0 = NULL, filter.end = NULL,
                     predict.end = NULL, minmax = c(-1.0e+30, 1.0e+30),
                     missed = NULL, np = NULL, plot = TRUE, ...)
{
# y                    # time series
  n <- length(y)
  ymean <- mean(y)

  ff <- f
  if (is.matrix(ff) == FALSE)
    ff <- as.matrix(ff)
  m <- dim(ff)[1]
  if (dim(ff)[2] != m)
    stop("Error : state transition matrix F(n) is invalid.")

  gg <- g
  if (is.matrix(gg) == FALSE)
    gg <- as.matrix(gg)
  if (dim(gg)[1] != m)
    stop("Error : matrix G(n) is invalid.")
  k <- dim(gg)[2]      # dimension of the system noise

  qq <- q
  if (is.matrix(qq) == FALSE)
    qq <- as.matrix(qq)
  if (dim(qq)[1] != k || dim(qq)[2] != k)
    stop("Error : matrix Q(n) is invalid.")

  if (is.null(x0))
    x0 <- rep(0.0e0, m)
  if (is.null(v0)) {
    v0 <- matrix(0.0e0, m, m )
    for (i in 1:m)
      v0[i, i] <- 10.0e-5
  }

  if (is.null(filter.end))
    filter.end <- n    # end point of filtering
  if (is.null(predict.end))
    predict.end <- n   # end point of prediction
  if (filter.end > predict.end)
    stop("'filter.end' is less than or equal to n (data length)")

  outmin <- minmax[1]   # lower limits of observations
  outmax <- minmax[2]   # upper limits of observations

  startp <- missed
  if (is.null(startp) || is.null(np)) {
    nmiss <- 0
    startp <- 0
    np <- 0
  } else {
    nmiss <- min(length(startp), length(np)) # number of missed intervals
    startp <- startp[1:nmiss]               # start position of missed intervals
    np <- np[1:nmiss]                       # number of missed observations
  }

  npe <- predict.end

  z <- .Fortran(C_smoothf,
             as.double(y),
             as.integer(n),
             as.integer(m),
             as.integer(k),
             as.double(ff),
             as.double(gg),
             as.double(h),
             as.double(qq),
             as.double(r),
             as.double(x0),
             as.double(v0),
             as.integer(filter.end),
             as.integer(npe),
             as.double(outmin),
             as.double(outmax),
             as.integer(nmiss),
             as.integer(startp),
             as.integer(np),
             xss = double(m * npe),
             vss = double(m * m * npe),
             lkhood = double(1),
             aic = double(1))

  xss <- array(z$xss, dim = c(m, npe)) + ymean
  vss <- array(z$vss, dim = c(m, m, npe))
  cov.smooth <- array(,dim = c(m, npe))
  for (i in 1:m)
    cov.smooth[i, ] <- vss[i, i, ]
  err <- NULL

  if (filter.end < n) {
    err <- rep(0, npe)
    nps <- filter.end + 1
    for (i in nps:npe)
      err[i] <- y[i] - xss[1, i]
  }

  if (nmiss != 0) {
    err <- rep(0, npe)
    for (i in 1:n)
      err[i] <- y[i] - xss[1, i]
  }

  tsmooth.out <- list(mean.smooth = xss, cov.smooth = cov.smooth, 
                      esterr = err, llkhood = z$lkhood, aic = z$aic)
  class(tsmooth.out) <- c("smooth")

  if (plot) {
    rdata <- deparse(substitute(y))
    eval(parse(text=paste(rdata, "<- y")))
    eval(parse(text=paste("plot.smooth(tsmooth.out,", rdata, ", ...)")))
  }

  return(tsmooth.out)
}

print.smooth <- function(x, ...) {
  message(gettextf("\n log-likelihood\t%12.3f", x$llkhood), domain = NA)
  message(gettextf(" aic\t\t%12.3f\n", x$aic), domain = NA)
}

plot.smooth <- function(x, rdata = NULL, ...)
{
  if (is.null(rdata) == FALSE)
    tsname <- deparse(substitute(rdata))

  old.par <- par(no.readonly = TRUE)
  new.mar <- old.par$mar
  new.mar[3] <- new.mar[3] + 1

  xss <- x$mean.smooth
  xss1 <- xss[1, ]
  cov <- x$cov.smooth
  npe <- dim(xss)[2]

  ts.atr <- tsp(rdata)
  if (is.null(ts.atr) == FALSE)
    xss1 <- ts(xss1, start = ts.atr[1], frequency = ts.atr[3])

  c1 <- xss1 + sqrt(cov[1, ])
  c2 <- xss1 - sqrt(cov[1, ])

  err <- x$esterr
  if (is.null(rdata) == FALSE && is.null(err)) {
    err <- rep(0, npe)
    n <- length(rdata)
    for (i in 1:n)
      err[i] <- rdata[i] - xss[1, i]
    if (min(err) == 0 && max(err) == 0)
      err <- NULL
  }

  if (is.null(err) == FALSE)
    par(mfcol=c(2,1))
  par(mar = new.mar, xpd=T)

  mtitle <- paste("Mean vectors of the smoother and standard deviation")
  ylim <- c(floor(min(rdata, c1, c2)), ceiling((max(rdata, c1, c2))))
  plot(c1, type='l', ylim = ylim, col=4, xlab="", ylab="", ...)
  par(new = TRUE)
  plot(c2, type='l', ylim = ylim, col=4, xlab="", ylab="", ...)
  par(new = TRUE)
  plot(xss1, type = 'l', ylim = ylim, col = 2, main = mtitle, xlab="", ylab="",
       ...)
  par(new = TRUE)

  pos.y <- par()$usr[4]
  if (is.null(rdata) == FALSE) {
    if (is.null(ts.atr) == TRUE) {
      xlim <- range(1, npe)
      plot(rdata, type = 'l', xlim = xlim, ylim = ylim, xlab = "", ylab = "", ...)
    } else {
      plot(rdata, type = 'l', ylim = ylim, xlab = "", ylab = "", ...)
    }
    lg3 <- legend("topright", legend = c("mean", "+/- SD", paste(tsname)),
                 lty = c(1, 1, 1), x.intersp = 0.3,ncol = 3, cex = 0.7, 
                 plot = FALSE)
    pos.y <- pos.y + lg3$rect$h
    legend(lg3$rect$left, pos.y, legend = c("mean", "+/- SD", paste(tsname)),
                 lty = c(1, 1, 1), col = c(2, 4, 1), x.intersp = 0.3, ncol = 3,
                 cex = 0.7, bty = "n")
  } else {
    lg2 <- legend("topright", legend = c("mean", "+/- SD"), lty = c(1, 1),
                 x.intersp = 0.3, ncol = 2, cex = 0.7, plot = FALSE)
    pos.y <- pos.y + lg2$rect$h
    legend(lg2$rect$left, pos.y, legend = c("mean", "+/- SD"), lty = c(1, 1),
           col = c(2, 4), x.intersp = 0.3,  ncol = 2, cex = 0.7, bty = "n")
  }

  if (is.null(err) == FALSE) {
    if (is.null(ts.atr) == FALSE)
      err <- ts(err, start = ts.atr[1], frequency = ts.atr[3])
    plot(err, type = 'h', main = "estimation error", xlab = "",
         ylab = "", ...)
  }
  par(old.par)
}
