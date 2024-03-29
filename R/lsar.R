# PROGRAM 8.1  LSAR1
lsar <- function (y, max.arorder = 20, ns0, plot = TRUE, ...)
{
  n <- length(y)         # Data length
# max.arorder            # Highest AR order
  lag <- max.arorder
# ns                     # Basic local span
  nb <- as.integer(n / ns0)
  nf0 <- 100             # Initial number of frequencies for computing spectrum

  z <- .Fortran(C_lsar1,
                as.double(y),
	            as.integer(n),
                as.integer(lag),
                as.integer(ns0),
                as.integer(nb),
                as.integer(nf0),
                ns = integer(nb),
                n0 = integer(nb),
                n1 = integer(nb),
                iflag = integer(nb),
                ms = integer(nb),
                sds = double(nb),
                aics = double(nb),
                mp = integer(nb),
                sdp = double(nb),
                aicp = double(nb),
                as = double(nb*lag),
                mfs = integer(nb),
                sig2s = double(nb),
                nnf = integer(nb),
                ier = integer(1))

  if (z$ier == 0) {
    model <- z$iflag
    sub <- list(start = z$n0, end = z$n1)
    as <- array(z$as, dim = c(lag, nb))
    nnf <- z$nnf

    spec <- sub.arsp(model, sub, as, z$mfs, z$sig2s, nnf)
    sdp <- z$sdp
    sdp[1] <- NA
    aicp <- z$aicp
    aicp[1] <- NA


    if (ns0 < max.arorder * 2)
      warning("the length of basic local span (ns0) is too small. 
ns0 should be at least twice as large as the maximum AR order")

    lsar.out <- list(model = model, ns = z$ns, span = sub, nf = nnf, ms = z$ms,
                     sds = z$sds, aics = z$aics, mp = z$mp, sdp = sdp, aicp = aicp,
                     spec = spec, tsname = deparse(substitute(y)))
    class(lsar.out) <- "lsar"

    if (plot)
      plot.lsar(lsar.out, ...)

    return(lsar.out)

  } else {
    stop("estimation error: the length of basic local span (ns0) is too small. 
ns0 should be at least twice as large as the maximum AR order")
  }
}

sub.arsp <- function(model, sub, as, mfs, sig2s, nf)
{
  st <- sub$start[1]
  ed <- NULL

  spec <- list()
  b <- 0
  l <- 0
  nb <- length(model)
  k <- 0
  for (i  in 1:(nb-1))
    if (model[i + 1] == 2) {
      a <- as[, i]
      m <- mfs[i]
      sig2 <- sig2s[i]
##      nnf <- nf[i]
      nnf <- max(nf[i], 100)

      z <- .Fortran(C_armasp,
                   as.double(a),
                   as.integer(m),
                   as.double(b),
                   as.integer(l),
                   as.double(sig2),
                   as.integer(nnf),
                   spec = double(nnf + 1))

      k <- k + 1
      spec[[k]] <- z$spec
      ed <- c(ed, sub$end[i])
      st <- c(st, sub$start[i + 1])
    }

  a <- as[, nb]
  m <- mfs[nb]
  sig2 <- sig2s[nb]
  nnf <- nf[nb]

  z <- .Fortran(C_armasp,
                as.double(a),
                as.integer(m),
                as.double(b),
                as.integer(l),
                as.double(sig2),
                as.integer(nnf),
                spec = double(nnf + 1))

  k <- k + 1
  spec[[k]] <- z$spec
  ed[k] <- sub$end[nb]

  spec <- list(subspec = spec, start = st, end = ed)
  return(spec = spec)
}

plot.lsar <- function(x, ...)
{
  spec <- x$spec$subspec
  st <- x$spec$start
  ed <- x$spec$end

  k <- length(spec)
  sp <- unlist(spec)
  ylim <- c(floor(min(sp)), ceiling(max(sp)))

  old.par <- par(no.readonly = TRUE)
  par(mfrow = c(3, 3), xaxs = "i", yaxs = "i")
  ic <- 0
  for (i in 1:k) {
#    if (ic == 4) {
    if (ic == 9) {
      par(ask = TRUE)
      ic <- 0
    }
    nf1 <- length(spec[[i]])
    xx <- rep(0, nf1)
    for (j in 1:nf1)
      xx[j] <- (j - 1) / (2 * (nf1-1))

    mtitle <- paste(st[i], " - ", ed[i])
    if (i == 1)
      mtitle <- paste(x$tsname, "\n", st[i], " - ", ed[i])
    plot(xx, spec[[i]], type = "l", main = mtitle, ylab = "", xlab = "f",
         ylim = ylim, ...)
    ic <- ic + 1
  }
  par(old.par)
}

print.lsar <- function(x, ...)
{
  nb <- length(x$ns)
  n0 <- x$span$start
  n1 <- x$span$end

  message(gettextf("<<< new data ( n = %d --- %d ) >>>", n0[1], n1[1]),
          domain = NA)
  message(gettextf(" initial model : NS = %5d", x$ns[1]), domain = NA)
  message(gettextf("\t\tms = %d   sds = %13.6e   aics = %12.3f",
          x$ms[1], x$sds[1], x$aics[1]), domain = NA)


  for (i in 2:nb) {
    np <- x$ns[i] + x$nf[i - 1]
    message(gettextf("\n<<< new data ( n = %d --- %d ) >>>", n0[i], n1[i]),
            domain = NA)
    message(gettextf(" switched model : ( nf = %d, ns = %d )", x$nf[i-1],
            x$ns[i]), domain = NA)
    message(gettextf("\t ms = %d   sds = %13.6e   aics = %12.3f", x$ms[i],
            x$sds[i], x$aics[i]), domain = NA)

    message(gettextf(" pooled model : ( np = %d )", np), domain = NA)
    message(gettextf("\t mp = %d   sdp = %13.6e   aicp = %12.3f", x$mp[i],
            x$sdp[i], x$aicp[i]), domain = NA)

    if (x$model[i] == 1)
      message("\n ***  pooled model accepted   ***")
    if (x$model[i] == 2)
      message("\n ***  switched model accepted   ***")
  }
}

# PROGRAM 8.2  LSAR2
lsar.chgpt <- function (y, max.arorder = 20, subinterval, candidate,
                        plot = TRUE, ...)
{
  n <- length(y)      # Data length
# max.arorder         # Highest AR order
  k <- max.arorder
  if (is.null(subinterval))
    subinterval <- c(1, n)
# n0 < n1 < n2 < ne <= n
  n0 <- subinterval[1]
  ne <- subinterval[2]      # [n0,ne] Time interval used for model fitting
  n1 <- candidate[1]
  n2 <- candidate[2]        # [n1,n2] Candidate for change point

  if (ne > n)
    stop("ne is less than or equal to n (data length).")
  n02k <- n0 + 2 * k
  if ((n02k == n1) || (n02k > n1))
    stop("\n the following condition 'n0 + 2k < n1' is not satisfied, where n0
    is the start point of 'subinterval', n1 the start point of 'candidate' and
    k is max.arorder.")
  n2k <- n2 + k
  if ((n2k == ne) || (n2k > ne))
    stop("\n the following condition 'n2 + k < ne' is not satisfied, where n2 is
    the end point of 'candidate', \n k is max.arorder and ne is the end point of
    'subinterval'.")

  m <- n2 - n1

  z <- .Fortran(C_lsar2,
                as.double(y),
                as.integer(n),
                as.integer(max.arorder),
                as.integer(n0),
                as.integer(n1),
                as.integer(n2),
                as.integer(ne),
                aics = double(m),
                aicmin = double(1),
                nmin = integer(1))

  change.pt <- n1 + z$nmin
  subint <- list(tsname = deparse(substitute(y)), st = n1+1, ed = n2,
                 y = y[(n1+1):n2])

  lsar.chgpt.out <- list(aic = z$aics, aicmin = z$aicmin,
                         change.point = change.pt, subint = subint )
  class(lsar.chgpt.out) <- "chgpt"

  if (plot) {
    plot.chgpt(lsar.chgpt.out, ...)
    invisible(lsar.chgpt.out)
  } else lsar.chgpt.out
}

plot.chgpt <- function(x, ...)
{
  cpt <- x$change.point
  n1 <- x$subint$st
  n2 <- x$subint$ed
  xx <- c(n1:n2)
  yy <- x$subint$y

  old.par <- par(no.readonly = TRUE)
  par(mfrow = c(2, 1), xaxs = "i")

  plot(xx, yy, type = "l", main = paste(x$subint$tsname, "sub-interval"),
       xlab = "", ylab = "y", cex.main = 1.0, ...)

  aicmin <- format(round(x$aicmin, 2), nsmall=2)
  mtitle <- paste("AICs of the model in the candidate range\n minimum aic =",
                  aicmin, " at ", cpt)
  plot(xx, x$aic, type = "l", main = mtitle, xlab = "time", ylab = "aic",
       cex.main = 1.0, ...)
  points(cpt, x$aicmin, col = 2, pch = 20)

  par(old.par)
}
