# PROGRAM 7.2 MARFIT
marfit <- function(y, lag = NULL)
{
  n <- dim(y)[1]              # data length
  l <- dim(y)[2]              # dimension of the observation
  if (is.null(lag))           # maximum lag of the cross-covariance function
    lag <- as.integer(2 * sqrt(n))

  z <- .Call("MarfitC",
             as.double(y),
             as.integer(n),
             as.integer(l),
             as.integer(lag))

  aic <- z[[3L]]
  mmin <- z[[4L]]
  amin <- array(z[[1L]], c(lag, l, l))
  vmin <- array(z[[2L]], c(l, l))
  arcoef <- array(, c(l, l, mmin))
  for (i in 1:mmin)
    arcoef[, , i] <- amin[i, , ]

  marfit.out <- list(maice.order = mmin, aic = aic, v = vmin, arcoef = arcoef)
  class(marfit.out) <- "maryule"
  return(marfit.out)
}

print.maryule <- function(x, ...)
{
  aic <- x$aic
  n <- length(aic)
  m <- dim(x$v)[1]
  morder <- x$maice.order
  aicmin <- aic[morder+1] 

  message("\n Order     AIC")
  for (i in 1:n) 
    message(gettextf("%5i  %12.3f", i-1, aic[i]), domain = NA)

  message(gettextf("\n Minimum AIC = %12.3f   attained at m = %5i", aicmin,
          morder), domain = NA)

  message("\n Innovation covariance matrix")
  for (i in 1:m) {
    for (j in 1:m)
      message(gettextf(" %13.6e", x$v[i,j]), appendLF=FALSE, domain = NA)
    message("")
  }
}

# PROGRAM 7.3 MARLSQ
marlsq <- function(y, lag = NULL)
{
  n <- dim(y)[1]             # data length
  l <- dim(y)[2]             # dimension of the observation
  if (is.null(lag))          # maximum lag of the cross-covariance function
    lag <- as.integer(2 * sqrt(n))

  z <- .Call("MarlsqC",
             as.double(y),
             as.integer(n),
             as.integer(l),
             as.integer(lag))

  arcoef <- array(z[[1L]], c(l, l, lag))
  v <- array(z[[2L]], c(l, l))
  lmax <- z[[3L]]
  aic <- z[[4L]]
  a <- arcoef[, , 1:lmax]

  marlsq.out <- list(maice.order = lmax, aic = aic, v = v, arcoef = a)
  class(marlsq.out) <- "marlsq"
  return(marlsq.out)
}

print.marlsq <- function(x, ...)
{

  message(gettextf("\n Total AIC = %12.3f", x$aic), domain = NA)
  message(gettextf(" order of the MAICE model = %5i", x$maice.order),
          domain = NA)

  m <- dim(x$v)[1]
  message("\n Innovation covariance matrix")
  for (i in 1:m) {
    for (j in 1:m)
      message(gettextf(" %13.6e", x$v[i,j]), appendLF=FALSE, domain = NA)
    message("")
  }
}
