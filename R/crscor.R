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

  y[is.na(y)] <- outmin

  z <- .Fortran(C_crscorf,
                as.double(y),
                as.integer(n),
                as.integer(id),
                as.integer(lag),
                as.double(outmin),
                as.double(outmax),
                cov = double(lag1 * id * id),
                cor = double(lag1 * id * id),
                mean = double(id))

  cov <- array(z$cov, c(lag1, id, id))
  cor <- array(z$cor, c(lag1, id, id))

  crscor.out <- list(cov = cov, cor = cor, mean = z$mean)
  class(crscor.out) <- "crscor"

  if (plot) {
    plot.crscor(crscor.out, colnames(y), ...)
    invisible(crscor.out)
  } else crscor.out
}

plot.crscor <- function(x, vnames, ...) {
  cor <- x$cor
  lag <- dim(cor)[1] - 1
  id <- dim(cor)[2]

  old.par <- par(no.readonly = TRUE)
    new.mar <- old.par$mar
    new.mar[3:4] <- new.mar[3:4] * 0.8
    new.mgp <- old.par$mgp
    new.mgp[1] <- new.mgp[1] * 0.8
    
    
  par(mfcol = c(id, id), xaxs = "i", yaxs = "i", oma = c(0, 1, 1, 0),
      mar = new.mar, mgp = new.mgp, cex.main = 0.9)

  for (j in 1:id) {
    for (i in 1:id) {
      xlab <- "";  ylab <- ""
      if (j == 1)
        ylab <- vnames[i]
      if (i == id)
        xlab <- "Lag"
      plot((0:lag), cor[, i, j], type = "l", ylim = c(-1, 1),
           xlab = xlab, ylab = ylab)
      if (i == 1)
        title(main = vnames[j], line = 1)

      par(new = TRUE)
      abline(h = 0)
    }
  }
  par(old.par)
}
