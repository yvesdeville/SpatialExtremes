rmaxstab <- function(n, coord, cov.mod = "gauss", grid = FALSE, ...){

  if (!(cov.mod %in% c("gauss","whitmat","cauchy","powexp","bessel")))
    stop("'cov.mod' must be one of 'gauss', 'whitmat', 'cauchy', 'powexp' or 'bessel'")

  dist.dim <- ncol(coord)
  
  if (is.null(dist.dim))
    dist.dim <- 1

  if (dist.dim > 2)
    stop("Currently this function is only available for R or R^2")

  if ((dist.dim == 1) && grid){
    warning("You cannot use 'grid = TRUE' in dimension 1. Ignored.")
    grid <- FALSE
  }

  if (cov.mod == "gauss"){
    if ((dist.dim == 1) && (!("var" %in% names(list(...)))))
      stop("You must specify 'var'")
    
    if ((dist.dim == 2) && (!all(c("cov11", "cov12", "cov22") %in% names(list(...)))))
      stop("You must specify 'cov11', 'cov12', 'cov22'")

    ##Get the model parameters
    var <- list(...)$var
    cov11 <- list(...)$cov11
    cov12 <- list(...)$cov12
    cov22 <- list(...)$cov22
    cov13 <- list(...)$cov13
    cov23 <- list(...)$cov23
    cov33 <- list(...)$cov33
  }

  else{
    if (!all(c("sill", "range", "smooth") %in% names(list(...))))
      stop("You must specify 'sill', 'range', 'smooth'")
    
    sill <- list(...)$sill
    range <- list(...)$range
    smooth <- list(...)$smooth
    alpha <- list(...)$alpha
    sigma2 <- list(...)$sigma2
  }

  if (dist.dim !=1){
    n.site <- nrow(coord)
    coord.range <- apply(coord, 2, range)
    center <- colMeans(coord.range)
    edge <- max(apply(coord.range, 2, diff))
  }

  else {
    n.site <- length(coord)
    coord.range <- range(coord)
    center <- mean(coord.range)
    edge <- diff(coord.range)
  }

  
  cov.mod <- switch(cov.mod, "gauss" = "gauss", "whitmat" = 1, "cauchy" = 2,"powexp" = 3,
                    "bessel" = 4)

  if (grid)
    ans <- double(n * n.site^dist.dim)

  else
    ans <- double(n * n.site)
  
  if (cov.mod == "gauss")
    ans <- switch(dist.dim,
                  .C("rsmith1d", as.double(coord), as.double(center), as.double(edge),
                     as.integer(n), as.integer(n.site), as.double(var), ans = ans,
                     PACKAGE = "SpatialExtremes")$ans,
                  .C("rsmith2d", as.double(coord), as.double(center), as.double(edge),
                     as.integer(n), as.integer(n.site), grid, as.double(cov11), as.double(cov12),
                     as.double(cov22), ans = ans, PACKAGE = "SpatialExtremes")$ans)
  else{
    if ((length(ans) / n) > 600)
      fun.name <- "rschlathertbm"

    else
      fun.name <- "rschlatherdirect"

    ans <- .C(fun.name, as.double(coord), as.integer(n), as.integer(n.site),
              as.integer(dist.dim), as.integer(cov.mod), grid, as.double(sill),
              as.double(range), as.double(smooth), ans = ans,
              PACKAGE = "SpatialExtremes")$ans
  }

  if (grid)
    ans <- array(ans, c(n.site, n.site, n))

  else
    ans <- matrix(ans, nrow = n, ncol = n.site)
  
  return(ans)
}
