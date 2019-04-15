# PROGRAM 10.1  ARMAFT
armafit <- function (y, ar.order, ar = NULL, ma.order, ma = NULL)
{
  n <- length(y)             # length of data
  m <- ar.order
  l <- ma.order
  mlmax <- max(m, l+1)
  if (is.null(ar) || is.null(ma)) {
    iparam <- 0
    ar <- rep(0, m)
    ma <- rep(0, l)
  } else {
    iparam <- 1
  }
    
  z <- .Call("armaft",
             as.double(y),
             as.integer(n),
             as.integer(m),
             as.integer(l),
             as.integer(mlmax),
             as.integer(iparam),
             as.double(ar),
             as.double(ma))

  ier <- z[[6L]]
  if (ier == 1)
    stop(" Matrix with zero row in decompose" )
  if (ier == 2)
    stop(" Singular matrix in decompose.zero divide in solve" )
  if (ier == 3)
    stop(" Convergence in impruv.matrix is nearly singular" )

  armafit.out <- list(sigma2 = z[[1L]], llkhood = z[[2L]], aic = z[[3L]],
                      arcoef = z[[4L]], macoef = z[[5L]])
  class(armafit.out) <- "armafit"

  return(armafit.out)
}

print.armafit <- function(x, ...)
{
  message(gettextf("\n sigma2\t\t%12.5e", x$sigma2), domain = NA)
  message(gettextf(" log-likelihood\t%12.3f", x$llkhood), domain = NA)
  message(gettextf(" aic\t\t%12.3f", x$aic), domain = NA)

 l <- length(x$arcoef)
  if (l > 0) {
    message("\n AR coefficients")
    for (i in 1:l)
      message(gettextf(" %12.3f\t", x$arcoef[i]), appendLF=FALSE, domain = NA)
    message("")
  }

  m <- length(x$macoef)
  if (m > 0) {
    message("\n MA coefficients")
    for (i in 1:m)
      message(gettextf(" %12.3f\t", x$macoef[i]), appendLF=FALSE, domain = NA)
    message("")
  }
}

