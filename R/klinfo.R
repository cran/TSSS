# PROGRAM 4.2 KLINFO
klinfo <- function(distg = 1, paramg = c(0, 1), distf = 1, paramf, xmax = 10)
{
  xmin <- -xmax

  z <- .Call("KlinfoC",
             as.integer(distg),
             as.double(paramg),
             as.integer(distf),
             as.double(paramf),
             as.double(xmin),
             as.double(xmax))

  klinfo.out <- list(nint = z[[1L]], dx = z[[2L]], KLI = z[[3L]],
                     gint = z[[4L]], xmax = xmax)
  class(klinfo.out) <- "klinfo"
  return(klinfo.out)
}

print.klinfo <- function(x, ...)
{
  n <- length(x$nint)
  message("\n  xmax\t    k\t    delta\t     KL info\t        gint")
  for (i in 1:n) 
    message(gettextf("%6.2f\t%5d\t%9.4f\t%12.8f\t%12.8f", x$xmax, x$nint[i],
                x$dx[i], x$KLI[i], x$gint[i]), domain = NA)
}
