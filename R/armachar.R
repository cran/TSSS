# PROGRAM 6.1
armachar <- function(arcoef = NULL, macoef = NULL, v, lag = 50, nf = 200,
                     plot = TRUE, ...)
{
  if (is.null(arcoef)) {        # AR coefficients
    arorder <- 0
    arcoef <- 0.0
  } else {
    arorder <- length(arcoef)   # AR order
  }
  if (is.null(macoef)) {        # MA coefficients
    maorder <- 0
    macoef <- 0.0
  } else {
    maorder <- length(macoef)   # MA order
  }
# v                             # innovation variance
# lag                           # maximum lag of autocovariance function
  kmax <- max(arorder, maorder, lag)
# nf                            # number of frequencies in evaluating spectrum

  z <- .Call("arma",
             as.integer(arorder),
             as.integer(maorder),
             as.double(arcoef),
             as.double(macoef),
             as.double(v),
             as.integer(lag),
             as.integer(kmax),    
             as.integer(nf))

  ier <- z[[7L]]
  if (ier == 1)
    stop(" Matrix with zero row in decompose" )
  if (ier == 2)
    stop(" Singular matrix in decompose.zero divide in solve" )
  if (ier == 3)
    stop(" Convergence in impruv.matrix is nearly singular" )

  impuls <- z[[1L]]
  acov <- z[[2L]]
  parcor <- z[[3L]]

  spec <- z[[4L]]
  if (arorder != 0)
    roota <- array(z[[5L]], dim = c(arorder, 2))
  if (maorder != 0)
    rootb <- array(z[[6L]], dim = c(maorder, 2))

  croot.ar <- NULL
  croot.ma <- NULL
  jer <- z[[8L]]
  if (jer == 1)
    warning(" AR : Non-convergence at polyrt\n")
  if (jer == 2)
    warning(" MA : Non-convergence at polyrt\n")
  if (jer == 3)
    warning(" Non-convergence at polyrt\n")

  if (arorder != 0) {
    croot.ar <- list()
    for (i in 1:arorder) {
      re <- roota[i, 1]
      im <- roota[i, 2]
      amp <- sqrt(re**2 + im**2)
      atan <- atan2(im, re)
      croot.ar[[i]] <- list(re = re, im = im, amp = amp, atan = atan,
                            degree = atan * 57.29577951)
    }
  }

  if (maorder != 0) {
    croot.ma <- list()
    for (i in 1:maorder) {
      re <- rootb[i, 1]
      im <- rootb[i, 2]
      amp <- sqrt(re**2 + im**2)
      atan <- atan2(im,re)
      croot.ma[[i]] <- list(re = re, im = im, amp = amp, atan = atan,
                            degree = atan*57.29577951)
    }
  }

  armachar.out <- list(impuls = impuls, acov = acov, parcor = parcor,
                      spec = spec, croot.ar = croot.ar, croot.ma = croot.ma)
  class(armachar.out) <- "arma"

  if (plot) {
    plot.arma(armachar.out, ...)
    invisible(armachar.out)
  } else armachar.out

}

plot.arma <- function(x, ...)
{
  impuls <- x$impuls
  acov <- x$acov
  parcor <- x$parcor
  spec <- x$spec
  croot.ar <- x$croot.ar
  croot.ma <- x$croot.ma

  old.par <- par(no.readonly = TRUE)
  par(xaxs = "i")
  if (is.null(croot.ar) == TRUE && is.null(croot.ma) == TRUE) {
    ier <- 3
    par(mfrow = c(2, 2), xaxs = "i")
  } else {
    ier <- 0
    par(mfrow = c(3, 2), xaxs = "i")
  }

  lag <- length(impuls) - 1
  x <- c(0:lag)
  ylim <- c(floor(min(impuls)), ceiling(max(impuls)))
  plot(x, impuls, type = "l", xlim = c(0, lag), ylim = ylim, xlab = "",
       ylab = "", ...)
  par(new = TRUE)
  plot(x, impuls, type = "h", xlim = c(0, lag), ylim = ylim, xlab = 'lag',
       ylab = 'impulse', ...)

  ylim <- c(floor(min(acov)), ceiling(max(acov)))
  plot(x, acov, type = "l", xlim = c(0, lag), ylim = ylim, xlab = "",
       ylab = "", ...)
  par(new = TRUE)
  plot(x, acov, type = "h", xlim = c(0, lag), ylim = ylim, xlab = 'lag',
       ylab = 'autocovariance', ...)

  plot(parcor, type = "l", xlim = c(0, lag), ylim = c(-1, 1), xlab = "",
       ylab = "", ...)
  par(new = TRUE)
  plot(parcor, type = "h", xlim = c(0, lag), ylim = c(-1, 1), xlab = 'lag',
       ylab = 'parcor', ...)

  k1 <- length(spec)
  k <- k1 - 1
  x <- rep(0, k1)
  for (i in 1:k1)
    x[i] <- (i - 1) / (2 * k)
  ylim <- c(floor(min(spec)), ceiling(max(spec)))
  plot(x, spec, type = "l", ylim = ylim, xlab = "frequency",
       ylab = "log spectrum", ...)

  if (ier == 0) {
    par(pty = "s")
    plot(x = c(-1.0, 1.0), y = c(0, 0), type = "l", xlim = c(-1.1, 1.1),
         ylim = c(-1.1, 1.1), xlab = "ARMA characteristic roots", ylab = "",
         axes = FALSE, ...)
    par(new = TRUE)
    plot(x = c(0, 0), y = c(-1.0, 1.0), type = "l", xlim = c(-1.1, 1.1),
         ylim = c(-1.1, 1.1), xlab = "", ylab = "", axes = FALSE, ...)
    symbols(x = 0, y = 0, circles = 1, xlim = c(-1.0, 1.0),
            ylim = c(-1.0, 1.0), inches = FALSE, add = TRUE, ...)

    if (is.null(croot.ar) == FALSE) {
      arorder <- length(croot.ar)
      for (i in 1:arorder) {     # characteristic roots of AR operator
        x1 <- croot.ar[[i]]$re
        y1 <- croot.ar[[i]]$im
        points(x1, y1, pch = 4, col = 6, cex = 0.8)
      }
    }

    if (is.null(croot.ma) == FALSE) {
      maorder <- length(croot.ma)
      for (i in 1:maorder) {     # characteristic roots of MA operator
        x2 <- croot.ma[[i]]$re
        y2 <- croot.ma[[i]]$im
        points(x2, y2, pch = 3, col = 4, cex = 0.8)
      }
    }

    par(xpd = T)
    legend(par()$usr[2], par()$usr[4], legend = c("AR", "MA"), col = c(6,4),
           pch = c(4, 3), cex = 0.8, bty="n")
  }

  par(old.par)
  invisible(x)
}