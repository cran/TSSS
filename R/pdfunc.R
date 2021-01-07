# PROGRAM 4.1
pdfunc <- function(model = "norm", mean=0, sigma2 = 1, mu = 0, tau2 = 1, shape,
                 lambda = 1, side = 1, df, xmin = 0, xmax = 1, plot = TRUE, ...)
{
  k <- 201

  if (model == "norm") {
    model.type <- 1
    param <- c(mean, sigma2, 0)

  } else if (model == "Cauchy") {
    model.type <- 2
    param <- c(mu, tau2, 0)

  } else if (model == "Pearson") {
    model.type <- 3
#    if ((shape < 0) || (shape == 0))
    if ((shape < 0.5) || (shape == 0.5))
      stop("'shape' is greater than 0.5")
    param <- c(mu, tau2, shape)

  } else if (model == "exp") {
    if (side == 1) {
      model.type <- 4
    } else if (side == 2) {
      model.type <- 0
    } else {
      stop("'side' is 1 or 2")
    }
    param <- c(lambda, side, 0)

  } else if (model == "Chi2") {
    model.type <- 5
    param <- c(df, 0, 0)

  } else if (model == "dexp") {
    model.type <- 6
    param <- rep(0, 3)

  } else if (model == "unif") {
    model.type <- 7
    param <- c(xmin, xmax, 0)

  } else {
    stop("the model type is invalid.")
  }

  z <- .Call("densty",
             as.integer(model.type),
             as.double(param),
             as.double(xmin),
             as.double(xmax),
             as.integer(k))

  pdfunc.out <- list(model = model, density = z[[1]], interval=c(xmin, xmax),
                     param=param)
  class(pdfunc.out) <- c("pdfunc")

  if (plot) {
    plot.pdfunc(pdfunc.out, ...)
    invisible(pdfunc.out)
  } else pdfunc.out
}

plot.pdfunc <- function(x, ...)
{
  k <- 201
  model <- x$model
  f <- x$density
  xmin <- x$interval[1]
  xmax <- x$interval[2]
  ddx <- xmax - xmin
  pa1 <- x$param[1]
  pa2 <- x$param[2]
  pa3 <- x$param[3]

  if (model == "norm") {
    mtitle = paste("Normal distribution N ( ", pa1, ",", pa2, ")")

  } else if (model == "Cauchy") {
    mtitle = paste("Cauchy distribution\n ( mu =", pa1, ", tau2 =", pa2, ")")

  } else if (model == "Pearson") {
    mtitle = paste("Pearson family of distributions\n ( mu =", pa1, ", tau2 =",
                    pa2, ", shape =", pa3, ")")
  } else if (model == "exp") {
    if (pa2 == 1)
     mtitle = paste("Exponential distribution\n ( lambda =", pa1, ")")
    if (pa2 == 2)
     mtitle = paste("Two-sided exponential distribution\n ( lambda =", pa1, ")")

  } else if (model == "Chi2") {
    mtitle = paste("Chi-square distribution\n ( k =", pa1, ")")

  } else if (model == "dexp") {
    mtitle = "Double exponential distribution"

  } else if (model == "unif") {
    mtitle = "Uniform distribution"
  }

  if (model != "unif") {
    xx <- xmin + c(0:(k-1)) * ddx/(k-1)
    plot(xx, f, type = "l", main = mtitle, xlab = "x", ylab = "density", ...)
  } else {
    xx <- xmin + c(1:(k-1)) * ddx/(k-1)
    plot(xx, f[2:k], type = "l", ylim=c(0, max(f)*2), main = mtitle,
         xlab = "x", ylab = "density", ...)
    lines(c(xx[1],xx[1]), c(0,f[2]))
    lines(c(xx[k-1],xx[k-1]), c(0,f[k]))
  }
}
