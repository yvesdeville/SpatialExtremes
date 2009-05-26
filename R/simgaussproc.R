rgp <- function(n, coord, cov.mod = "powexp", mean = 0, nugget = 0, 
                sill = 1, range = 1, smooth = 1, grid = FALSE, ...){

  if (grid && is.null(dim(coord)))
    stop("'grid' cannot be 'TRUE' if you specify univariate coordinates")

  if (!(cov.mod %in% c("whitmat","cauchy","powexp","bessel")))
    stop("''cov.mod'' must be one of 'whitmat', 'cauchy', 'powexp', 'bessel'")

  if (cov.mod == "whitmat")
    cov.mod.num <- 1
  if (cov.mod == "cauchy")
    cov.mod.num <- 2
  if (cov.mod == "powexp")
    cov.mod.num <- 3
  if (cov.mod == "bessel")
    cov.mod.num <- 4
  
  ##Identify the most accurate method for simulation
  if (grid && (nrow(coord)^ncol(coord) > 500))
    method <- "tbm"

  else if (nrow(coord) > 500)
    method <- "tbm"

  else
    method <- "exact"

  gp <- switch(method,
               "tbm" = .tbmgp(n, coord, cov.mod.num, nugget, sill, range,
                 smooth, grid, ...),
               "exact" = .exactgp(n, coord, cov.mod.num, nugget, sill, range,
                 smooth, grid, ...))

  return(mean + gp)
}

.exactgp <- function(n, coord, cov.mod, nugget, sill, range, smooth, grid){
  dist.dim <- ncol(coord)
  n.site <- nrow(coord)

  if (is.null(dist.dim)){
    dist.dim <- 1
    n.site <- length(coord)
  }

  if (grid)
    n.effsite <- n.site^dist.dim

  else
    n.effsite <- n.site

  gp <- .C("direct", as.integer(n), as.integer(n.site), grid, as.integer(cov.mod), 
           as.double(coord), as.integer(dist.dim), as.double(nugget), 
           as.double(sill), as.double(range), as.double(smooth), ans = double(n.effsite * n), 
           PACKAGE = "SpatialExtremes")$ans
  
  if (grid){
    if (dist.dim == 2)
      gp <- array(gp, c(n.site, n.site, n))

    else
      gp <- array(gp, c(n.site, n.site, n.site, n))
  }

  else
    gp <- matrix(gp, ncol = n.site, nrow = n)
  
  return(gp)
}

.tbmgp <- function(n, coord, cov.mod, nugget, sill, range, smooth, grid,
                   nlines = 1000){

  n.site <- nrow(coord)
  dim <- ncol(coord)

  if (grid)
    ans <- double(n.site^dim * n)

  else
    ans <- double(n.site * n)
  
  gp <- .C("tbm", as.integer(n), as.integer(n.site), as.integer(dim),
           as.integer(cov.mod), grid, as.double(coord), as.double(nugget),
           as.double(sill), as.double(range), as.double(smooth),
           as.integer(nlines), ans = ans, PACKAGE = "SpatialExtremes")$ans

  if (grid){
    if (dim == 2)
      gp <- array(gp, c(n.site, n.site, n))

    else
      gp <- array(gp, c(n.site, n.site, n.site, n))
  }

  else
    gp <- matrix(gp, nrow = n, ncol = n.site)

  return(gp)
}
