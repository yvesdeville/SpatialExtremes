rmaxstab <- function(n, coord, cov.mod = "gauss", grid = FALSE,
                     control = list(), ...){

  if (!(cov.mod %in% c("gauss","whitmat","cauchy","powexp","bessel",
                       "iwhitmat", "icauchy", "ipowexp", "ibessel",
                       "gwhitmat", "gcauchy", "gpowexp", "gbessel")))
    stop("'cov.mod' must be one of 'gauss', '(i/g)whitmat', '(i/g)cauchy', '(i/g)powexp' or '(i/g)bessel'")

  if (!is.null(control$method) && !(control$method %in% c("exact", "tbm")))
    stop("the argument 'method' for 'control' must be one of 'exact' and 'tbm'")

  if (cov.mod == "gauss")
    model <- "Smith"

  else if (cov.mod %in% c("whitmat","cauchy","powexp","bessel"))
    model <- "Schlather"

  else {
    cov.mod <- substr(cov.mod, 2, 10)

    if (substr(cov.mod, 1, 1) == "i")
      model <- "iSchlather"

    else
      model <- "Geometric"
  }

  dist.dim <- ncol(coord)
  
  if (is.null(dist.dim))
    dist.dim <- 1

  if (dist.dim > 2)
    stop("Currently this function is only available for R or R^2")

  if ((dist.dim == 1) && grid){
    warning("You cannot use 'grid = TRUE' in dimension 1. Ignored.")
    grid <- FALSE
  }

  if (model == "Smith"){
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

  else if (model == "Schlather"){
    if (!all(c("sill", "range", "smooth") %in% names(list(...))))
      stop("You must specify 'sill', 'range', 'smooth'")
    
    sill <- list(...)$sill
    range <- list(...)$range
    smooth <- list(...)$smooth
  }

  else if (model == "iSchlather"){
    if (!all(c("alpha", "sill", "range", "smooth") %in% names(list(...))))
      stop("You must specify 'alpha', 'sill', 'range', 'smooth'")
    
    sill <- list(...)$sill
    range <- list(...)$range
    smooth <- list(...)$smooth
    alpha <- list(...)$alpha
  }

  else {
    if (!all(c("sigma2", "sill", "range", "smooth") %in% names(list(...))))
      stop("You must specify 'sigma2', 'sill', 'range', 'smooth'")
    
    sill <- list(...)$sill
    range <- list(...)$range
    smooth <- list(...)$smooth
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

  
  cov.mod <- switch(cov.mod, "gauss" = "gauss", "whitmat" = 1, "cauchy" = 2,
                    "powexp" = 3, "bessel" = 4)

  if (grid)
    ans <- double(n * n.site^dist.dim)

  else
    ans <- double(n * n.site)
  
  if (model == "Smith")
    ans <- switch(dist.dim,
                  .C("rsmith1d", as.double(coord), as.double(center), as.double(edge),
                     as.integer(n), as.integer(n.site), as.double(var), ans = ans,
                     PACKAGE = "SpatialExtremes")$ans,
                  .C("rsmith2d", as.double(coord), as.double(center), as.double(edge),
                     as.integer(n), as.integer(n.site), grid, as.double(cov11), as.double(cov12),
                     as.double(cov22), ans = ans, PACKAGE = "SpatialExtremes")$ans)
  
  else if (model == "Schlather"){

    if (is.null(control$method)){
      if ((length(ans) / n) > 600)
        method <- "tbm"
        
      else
        method <- "direct"
    }

    else
      method <- control$method

    if (is.null(control$uBound))
      uBound <- 3.5

    else
      uBound <- control$uBound

    if (method == "direct")
      ans <- .C("rschlatherdirect", as.double(coord), as.integer(n), as.integer(n.site),
                as.integer(dist.dim), as.integer(cov.mod), grid, as.double(sill),
                as.double(range), as.double(smooth), as.double(uBound), ans = ans,
                PACKAGE = "SpatialExtremes")$ans

    else {
      if (is.null(control$nlines))
        nlines <- 1000

      else
        nlines <- control$nlines
      
      ans <- .C("rschlathertbm", as.double(coord), as.integer(n), as.integer(n.site),
                as.integer(dist.dim), as.integer(cov.mod), grid, as.double(sill),
                as.double(range), as.double(smooth), as.double(uBound), as.integer(nlines),
                ans = ans, PACKAGE = "SpatialExtremes")$ans
    }
  }

  else if (model == "Geometric"){

    if (is.null(control$method)){
      if ((length(ans) / n) > 600)
        method <- "tbm"
    
      else
        method <- "direct"
    }

    else
      method <- control$method

    if (is.null(control$uBound))
      uBound <- exp(3.5 * sqrt(sigma2) - 0.5 * sigma2)

    else
      uBound <- control$uBound

    if (method == "direct")
      ans <- .C("rgeomdirect", as.double(coord), as.integer(n), as.integer(n.site),
                as.integer(dist.dim), as.integer(cov.mod), grid, as.double(sigma2),
                as.double(sill), as.double(range), as.double(smooth),
                as.double(uBound), ans = ans, PACKAGE = "SpatialExtremes")$ans

    else {

      if (is.null(control$nlines))
        nlines <- 1000

      else
        nlines <- control$nlines
      
      ans <- .C("rgeomtbm", as.double(coord), as.integer(n), as.integer(n.site),
                as.integer(dist.dim), as.integer(cov.mod), grid, as.double(sigma2),
                as.double(sill), as.double(range), as.double(smooth), as.double(uBound),
                as.integer(nlines), ans = ans, PACKAGE = "SpatialExtremes")$ans
    }
  }

  else
    stop("not implemented yet")

  if (grid)
    ans <- array(ans, c(n.site, n.site, n))

  else
    ans <- matrix(ans, nrow = n, ncol = n.site)
  
  return(ans)
}
