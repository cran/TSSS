# PROGRAM 4.2 KLINFO
klinfo <- function(distg = 1, paramg = c(0, 1), distf = 1, paramf, xmax = 10)
{
  xmin <- -xmax

  z <- .Fortran(C_klinfof,
                as.integer(distg),
                as.double(paramg),
                as.integer(distf),
                as.double(paramf),
                as.double(xmin),
                as.double(xmax),
                nint = integer(4),
                dx = double(4),
                fkli = double(4),
                gint = double(4))
			 
  klinfo.out <- list(nint = z$nint, dx = z$dx, KLI = z$fkli,
                     gint = z$gint, xmax = xmax)
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
