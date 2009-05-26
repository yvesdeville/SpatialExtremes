covariance <- function(fitted, sill, range, smooth, cov.mod = "whitmat",
                       plot = TRUE, dist, xlab, ylab, ...){

  if (!missing(fitted)){
    cov.mod <- fitted$cov.mod
    smooth <- fitted$param["smooth"]
    range <- fitted$param["range"]
    sill <- fitted$param["sill"]
  }

  if (cov.mod == "gauss")
    stop("''covariance'' is not implemented for the Smith's model")

  if (!(cov.mod %in% c("whitmat", "cauchy", "powexp", "bessel")))
    stop("Invalid covariance model. ''cov.mod'' must be one of 'whitmat', 'cauchy', 'powexp', 'bessel'")
  
  if (cov.mod == "whitmat"){
    if ((smooth <= 0) || (range <= 0) || (smooth > 150) || (sill <= 0))
      stop("invalid parameter for the whittle-matern covariance function")
    
    cov.fun <- function(dist) {
      idx <- dist == 0
      ans <- rep(sill, length(dist))
      ans[!idx] <- sill * 2^(1-smooth) / gamma(smooth) * (dist[!idx] / range)^smooth *
        besselK(dist[!idx] / range, smooth)
      return(ans)
    }
  }

  if (cov.mod == "cauchy"){
    if ((smooth <= 0) || (range <= 0) || (sill <= 0))
      stop("invalid parameter for the cauchy covariance function")
    
    cov.fun <- function(dist) sill * (1 + (dist / range)^2)^-smooth
  }

  if (cov.mod == "powexp"){
    if ((smooth < 0) || (smooth > 2) || (range <= 0) || (sill <= 0))
      stop("invalid parameter for the powered exponential covariance function")

    cov.fun <- function(dist) sill * exp(-(dist / range)^smooth)
  }

  if (cov.mod == "bessel"){
    if ((range <= 0) || (sill <= 0))
      stop("invalid parameter for the Bessel covariance function")

    cov.fun <- function(dist) sill * (2 * range / dist)^smooth * gamma(smooth + 1) *
      besselJ(dist / range, smooth)
  }

  if (plot){

    if (missing(xlab))
      xlab <- "h"

    if (missing(ylab))
      ylab <- expression(rho(h))

    if (is.null(list(...)$xlim)){
      tmp.fun <- function(dist) (cov.fun(dist) - 0.05)^2
      xlimsup <- optimize(tmp.fun, c(1e-6, 10 * sqrt(sill) * range))$minimum
    }

    else
      xlimsup <- list(...)$xlim[2]
    
    plot(cov.fun, from = 0, to = xlimsup, xlab = xlab, ylab = ylab, ...)
  }

  if (!missing(dist)){
    cov.val <- cov.fun(dist)
    return(list(cov.fun = cov.fun, cov.val = cov.val))
  }

  else
    invisible(cov.fun)  
}
