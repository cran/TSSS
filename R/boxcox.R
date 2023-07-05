# PROGRAM 4.4 BOXCOX
boxcox <- function(y, plot = TRUE, ...)
{
  n <- length(y)             # length of data
  for (i in 1:n)
    if ((y[i] == 0) || (y[i] < 0))
      stop("Log-transformation cannot be applied to zeros and nagative numbers")

  z <- .Fortran(C_boxcoxf,
                as.double(y),
                as.integer(n),
                aiczt = double(21),
                ffzt = double(21),
                aicz = double(21),
                ffz = double(21),
                mean = double(21),
                var = double(21),
                zz = double(n))

  boxcox.out <- list(y, mean = z$mean, var = z$var, aic = z$aicz,
                     llkhood = z$ffz, z = z$zz, aic.z = z$aiczt,
                     llkhood.z = z$ffzt)

  class(boxcox.out) <- "boxcox"

  if (plot) {
    rdata <- deparse(substitute(y))
    eval(parse(text=paste(rdata, "<- y")))
    eval(parse(text=paste("plot.boxcox(boxcox.out,", rdata, ", ...)")))
  }

  return(boxcox.out)
}

print.boxcox <- function(x, ...)
{
  message("\n lambda    aic'         LL'         aic          LL        mean",
          appendLF=FALSE)
  message("           variance")

  aicm <- x$aic.z[1]
  am <- 1
  for (i in 1:21) {
    a <- (-i + 11) / 10
    message(gettextf("%6.2f %11.2f %11.2f %11.2f %11.2f   %e   %e", a,
            x$aic.z[i], x$llkhood.z[i], x$aic[i], x$llkhood[i], x$mean[i],
            x$var[i]), domain = NA)
    if (x$aic.z[i] < aicm) {
      aicm <- x$aic.z[i]
      am <- a
    } 
  }
  message(gettextf("\n lambda = %6.2f\t AIC' minimum = %11.2f\n", am, aicm),
          domain = NA)
}

plot.boxcox <- function(x, rdata = NULL, ...)
{
  ts.atr <- tsp(rdata)
  old.par <- par(no.readonly = TRUE)
  par(xaxs = "i")

  if (is.null(rdata) == FALSE) {
    tsname <- deparse(substitute(rdata))
    mtitle <- paste(tsname)
    par(mfcol = c(2, 1))
    plot(rdata, type = "l", main = mtitle, xlab = "", ylab = "y", ...)
  }

  tdata <- x$z
  if (is.null(ts.atr) == FALSE)
    tdata <- ts(tdata, start = ts.atr[1], frequency = ts.atr[3])
  plot(tdata, type = "l", main = "Transformed data", xlab = "",
       ylab = "log y", ...)
  par(old.par)
}
