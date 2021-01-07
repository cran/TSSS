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

  m <- length(x$arcoef)
  if (m > 0) {
    message("\n AR coefficients")
    for (i in 1:m)
      message(gettextf(" %12.3f\t", x$arcoef[i]), appendLF=FALSE, domain = NA)
    message("")
  }

  l <- length(x$macoef)
  if (l > 0) {
    message("\n MA coefficients")
    for (i in 1:l)
      message(gettextf(" %12.3f\t", x$macoef[i]), appendLF=FALSE, domain = NA)
    message("")
  }
}

#------------------------------------------------------------------------------

armafit2 <- function (y, ar.order, ma.order)
{
  n <- length(y)             # length of data
  mmax <- ar.order
  lmax <- ma.order
  mlmax <- max(mmax, lmax, 1) + 1

  z <- .Call("armaft2",
             as.double(y),
             as.integer(n),
             as.integer(mmax),
             as.integer(lmax),
             as.integer(mlmax))

  ier <- z[[6L]]
  if (ier == 1)
    stop(" Matrix with zero row in decompose" )
  if (ier == 2)
    stop(" Singular matrix in decompose.zero divide in solve" )
  if (ier == 3)
    stop(" Convergence in impruv.matrix is nearly singular" )

  sig2 <- array(z[[1L]], c(mmax + 1, lmax + 1))
  flk <- array(z[[2L]], c(mmax + 1, lmax + 1))
  aic <- array(z[[3L]], c(mmax + 1, lmax + 1))
  ar <- array(z[[4L]], c(mmax, mmax + 1, lmax + 1))
  ma <- array(z[[5L]], c(lmax, mmax + 1, lmax + 1))
  ar <- aperm(ar, c(2, 3, 1))
  ma <- aperm(ma, c(2, 3, 1))

  mm1 <- mmax + 1
  ll1 <- lmax + 1
  coef <- list()
  length(coef) <- mm1
  for (i in 1:mm1) {
    coef[[i]] <- list()
    length(coef[[i]]) <- ll1
    for (j in 1:ll1) {
      if (i == 1) {
        arcoef <- NULL
      } else {
        arcoef <- ar[i, j, 1:(i-1)]
      }   
      if (j == 1) {
        macoef <- NULL
      } else {
        macoef <- ma[i, j, 1:(j-1)]
      }
      coef[[i]][[j]] <- list(ar = arcoef, ma = macoef)
    }
 }

  aic.min <- min(aic, na.rm=T)
  morder <- which(aic == aic.min, arr.ind=TRUE)
  m <- morder[1]
  l <- morder[2]
  m1 <- m - 1
  l1 <- l - 1

  armafit2.out <- list(aicmin = aic.min, maice.order = c(m1, l1), sigma2 = sig2,
                     llkhood = flk, aic = aic, coef = coef)

  class(armafit2.out) <- "armafit2"
  return(armafit2.out)
}

print.armafit2 <- function(x, ...)
{
  mmax <- dim(x$sigma2)[1] - 1
  lmax <- dim(x$sigma2)[2] - 1

  sigma2 <- x$sigma2
  llkhood <- x$llkhood
  aic <- x$aic
  rownames(sigma2) <- c(0:mmax)
  colnames(sigma2) <- c(0:lmax)
  rownames(llkhood) <- c(0:mmax)
  colnames(llkhood) <- c(0:lmax)
  rownames(aic) <- c(0:mmax)
  colnames(aic) <- c(0:lmax)

  m <- x$maice.order[1]
  l <- x$maice.order[2]
  m1 <- m + 1
  l1 <- l + 1

  message(gettextf("\n minimum AIC %12.3f at order( %d, %d )", x$aicmin, m, l),
          domain = NA)
  if (m > 0) {
    message(" AR coefficients")
    arcoef <- x$coef[[m1]][[l1]]$ar
    for (i in 1:m) {
      message(gettextf(" %12.3f\t", arcoef[i]), appendLF=FALSE, domain = NA)
      if ((i %% 5 == 0) || i == m)
        message("")
    }
  }
  if (l > 0) {
    message(" MA coefficients")
    macoef <- x$coef[[m1]][[l1]]$ma
    for (i in 1:l) {
      message(gettextf(" %12.3f\t", macoef[i]), appendLF=FALSE, domain = NA)
      if ((i %% 5 == 0) || i == l)
        message("")
    }
  }
  message(gettextf("\n sigma2(m,l) (0 <= m <= %d, 0 <= l <= %d)", mmax, lmax), domain = NA)
  message(paste0(capture.output(sigma2), collapse = "\n"))
  message(gettextf("\n llkhood(m,l) (0 <= m <= %d, 0 <= l <= %d)", mmax, lmax), domain = NA)
  message(paste0(capture.output(llkhood), collapse = "\n"))
  message(gettextf("\n aic(m,l) (0 <= m <= %d, 0 <= l <= %d)", mmax, lmax), domain = NA)
  message(paste0(capture.output(aic), collapse = "\n"))
  message("")
}

