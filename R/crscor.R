
# PROGRAM 2.1
crscor <- function (y, lag = NULL, outmin = NULL, outmax = NULL, plot = TRUE, ...)
{
  n <- dim(y)[1]
  id <- dim(y)[2]
  if (is.null(lag))
    lag <- as.integer(2 * sqrt(n))    # maximum lag
  lag1 <- lag + 1

  if (is.null(outmin)) {
    outmin <- rep(-1.0e+30, id)
  } else {
    if (length(outmin) != id)
      stop("specify the values of 'outmin' for each dimension")
  }
  if (is.null(outmax)) {
    outmax <- rep(1.0e+30, id)
  } else {
    if (length(outmax) != id)
      stop ("specify the values of 'outmax' for each dimension")
  }

  z <- .Call("CrscorC",
             as.double(y),
             as.integer(n),
             as.integer(id),
             as.integer(lag),
             as.double(outmin),
             as.double(outmax))

  cov <- array(z[[1]], c(lag1, id, id))
  cor <- array(z[[2]], c(lag1, id, id))
  ymean <- z[[3L]]

  crscor.out <- list(cov = cov, cor = cor, mean = ymean)
  class(crscor.out) <- "crscor"

  if (plot) {
    plot.crscor(crscor.out, ...)
    invisible(crscor.out)
  } else crscor.out
}

plot.crscor <- function(x, ...) {
  cor <- x$cor
  lag <- dim(cor)[1] - 1
  id <- dim(cor)[2]

  old.par <- par(no.readonly = TRUE)
  par(mfcol = c(id, id), xaxs = "i", yaxs = "i")
  for (j in 1:id)
    for (i in 1:id) {
      plot((0:lag), cor[, i, j], type = "l", ylim = c(-1, 1),
           ylab = paste("cor [", i, ",", j, "]"), xlab = "Lag", ...)
      par(new = TRUE)
      abline(h = 0)
    }
  par(old.par)
}
