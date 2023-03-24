# PROGRAM 6.2 MARSPC
marspc <- function(arcoef, v, plot = TRUE, ...)
{
  m <- dim(arcoef)[3]      # ar order
  l <- dim(v)[1]           # dimension of the time series

  nf <- 200                # number of frequencies
  nf1 <- nf + 1
  llnf1 <- l * l * nf1

  z <- .Fortran(C_marspcf,
             as.integer(m),
             as.integer(l),
             as.double(arcoef),
             as.double(v),
             as.integer(nf),
             p = complex(llnf1),
             amp = double(llnf1),
             ang = double(llnf1),
             coh = double(llnf1),
             fnc = double(llnf1),
             frnc = double(llnf1))

  p <- array(z$p, dim = c(nf1, l, l))
  amp <- array(z$amp, dim = c(nf1, l, l))
  ang <- array(z$ang, dim = c(nf1, l, l))
  coh <- array(z$coh, dim = c(nf1, l, l))
  fnc <- array(z$fnc, dim = c(nf1, l, l))
  frnc <- array(z$frnc, dim = c(nf1, l, l))

  marspc.out <- list(spec = p, amp = amp, phase = ang, coh = coh, power = fnc,
                     rpower = frnc)

  if (plot) {
    plot.marspc(marspc.out, ...)
    invisible(marspc.out)
  } else marspc.out
}

plot.marspc <- function(x, ...)
{
  spec <- x$spec
  amp <- x$amp
  phase <- x$phase
  coh <- x$coh
  fnc <- x$power
#  frnc <- x$rpower

  old.par <- par(no.readonly = TRUE)

  d <- dim(spec)[3]
  d1 <- d - 1
  nf1 <- dim(spec)[1]
  fsp <- array(0, dim = c(nf1, d, d))

  x <- rep(0, nf1)
  for (i in 1:nf1)
    x[i] <- (i - 1) / (2 * (nf1 - 1))

# cross spectra

  for (i in 1:d)
    for(k in 1:nf1) 
      fsp[k, i, i] <- log10(Re(spec[k, i, i]))

  for (i in 1:d1) 
    for (j in (i+1):d)
      for(k in 1:nf1) {
        fsp[k, i, j] <- log10(amp[k, i, j])
        fsp[k, j, i] <- phase[k, i, j]
    }

  min1 <- min(fsp[, 1, 1])
  max1 <- max(fsp[, 1, 1])
  min2 <- min(fsp[, 1, 2])
  max2 <- max(fsp[, 1, 2])
  min3 <- min(fsp[, 2, 1])
  max3 <- max(fsp[, 2, 1])

  for (i in 2:d) {
    min1 <- min(min1, fsp[, i, i])
    max1 <- max(max1, fsp[, i, i])
  }
  if (d > 2) {
    for (i in 1:d1) 
      for (j in (i+1):d) {
          min2 <- min(min2, fsp[, i, j])
          max2 <- max(max2, fsp[, i, j])
          min3 <- min(min3, fsp[, j, i])
          max3 <- max(max3, fsp[, j, i])
      }
  }

  ylim1 <- c(floor(min1), ceiling(max1))
  ylim2 <- c(floor(min2), ceiling(max2))
  ylim3 <- c(floor(min3), ceiling(max3))

  par(mfrow = c(d, d), xaxs = "i", yaxs = "i")
  for (i in 1:d)
    for (j in 1:d) 
      if (i == j) {
        plot(x, fsp[, i, i], type = "l", ylim = ylim1, 
             xlab = "Frequency", ylab = paste("Spec (", i, ",", i, ")"), ...)
      } else if (i < j) {
        plot(x, fsp[, i, j], type = "l", ylim = ylim2,
             xlab = "Frequency", ylab = paste("Amp (", i, ",", j, ")"), ...)
      } else {
        plot(x, fsp[, i, j], type = "l", ylim = ylim3,
             xlab = "Frequency", ylab = paste("Phase (", i, ",", j, ")"), ...)
      }

# Power spectre and coherency
  par(ask = TRUE)
  par(mfrow = c(d, d), xaxs = "i", yaxs = "i")
  ymax <- max(coh)
  for (i in 1:d)
    for (j in i:d) {
      if (j != 1) par(mfg = c(i, j, d, d))
      if(i == j) {
        plot(x, fsp[, i, i], type = "l", ylim = ylim1,
             xlab = "Frequency", ylab = paste("Spec (", i, ",", i, ")"), ...)
      } else {
        plot(x, coh[, i, j], type = "l", ylim = c(0, ymax),
             xlab = "Frequency", ylab = paste("Coh (", i, ",", j, ")"), ...)
      }
    }

  par(ask=TRUE)
  par(mfrow = c(d,2))
  ymin <- 0
  for(i in 1:d) {
    ymax <- max(fnc[,i,])*1.1
    for(j in 1:d) {
      if( i == 1 )
        plot(x, fnc[,i,j], type='l', ylim=c(ymin,ymax), xlab="", ylab="",
             main="Power contribution", ...)
      if( i != 1 )
        plot(x, fnc[,i,j], type='l', ylim=c(ymin,ymax), xlab="", ylab="", ...)
      par(new=TRUE)
    }
    par(new=FALSE)

    for(j in 1:d) {
      rnc <- fnc[, i, j] / fnc[, i, d]
      if (i == 1)
        plot(x, rnc, type = "l", ylim = c(0, 1), xlab = "", ylab = "",
             main = "Relative power contribution", ...)
      if (i != 1)
        plot(x, rnc, type = "l", ylim = c(0, 1), xlab = "", ylab = "", ...)
      par(new = TRUE)
    }
    par(new = FALSE)
  }
  par(old.par)
}
