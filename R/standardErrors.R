.smithstderr <- function(par, data, distVec, loc.dsgn.mat,
                         scale.dsgn.mat, shape.dsgn.mat, fit.marge,
                         std.err.type = "score", fixed.param, param.names,
                         iso = TRUE, weights){

  ##data is a matrix with each column corresponds to one location
  ##distVec is the a matrix giving the "distance vector" for each pair
  ##(1 row = 1 station)
  n.site <- ncol(data)
  n.obs <- nrow(data)
  n.pairs <- n.site * (n.site - 1) / 2
  dist.dim <- ncol(distVec)
  n.param <- length(param.names)

  if (iso){
    if (dist.dim == 2)
      n.param <- n.param + 2

    else
      n.param <- n.param + 5
  }
  
  cov11 <- par["cov11"]
  cov12 <- par["cov12"]
  cov22 <- par["cov22"]
  
  if (dist.dim == 3){
    cov13 <- par["cov13"]
    cov23 <- par["cov23"]
    cov33 <- par["cov33"]
  }
        
  if (fit.marge){

    n.loccoeff <- ncol(loc.dsgn.mat)
    n.scalecoeff <- ncol(scale.dsgn.mat)
    n.shapecoeff <- ncol(shape.dsgn.mat)

    loc.idx <- which(substr(names(par), 1, 3) == "loc")
    scale.idx <- which(substr(names(par), 1, 5) == "scale")
    shape.idx <- which(substr(names(par), 1, 5) == "shape")

    loc.param <- par[loc.idx]
    scale.param <- par[scale.idx]
    shape.param <- par[shape.idx]
  }

  else {
    n.loccoeff <- 1
    n.scalecoeff <- 1
    n.shapecoeff <- 1
  
    loc.param <- 1
    scale.param <- 1
    shape.param <- 1
  }

  if (is.null(weights)){
    if (dist.dim == 2)
      std.err <- .C("smithstderr", as.double(data), as.double(distVec), as.integer(n.site),
                    as.integer(n.obs), as.double(loc.dsgn.mat), as.integer(n.loccoeff),
                    as.double(scale.dsgn.mat), as.integer(n.scalecoeff), as.double(shape.dsgn.mat),
                    as.integer(n.shapecoeff), as.double(loc.param), as.double(scale.param),
                    as.double(shape.param), as.double(cov11), as.double(cov12),
                    as.double(cov22), fit.marge, hess = double(n.obs * n.param * n.pairs),
                    grad = double(n.obs * n.param), PACKAGE = "SpatialExtremes")
    
    else
      std.err <- .C("smithgrad3d", as.double(data), as.double(distVec), as.integer(n.site),
                    as.integer(n.obs), as.double(loc.dsgn.mat), as.integer(n.loccoeff),
                    as.double(scale.dsgn.mat), as.integer(n.scalecoeff), as.double(shape.dsgn.mat),
                    as.integer(n.shapecoeff), as.double(loc.param), as.double(scale.param),
                    as.double(shape.param), as.double(cov11), as.double(cov12), as.double(cov13),
                    as.double(cov22), as.double(cov23), as.double(cov33), fit.marge,
                    hess = double(n.obs * n.param * n.pairs), grad = double(n.obs * n.param),
                    PACKAGE = "SpatialExtremes")
  }

  else{
   if (dist.dim == 2)
      std.err <- .C("wsmithstderr", as.double(data), as.double(distVec), as.integer(n.site),
                    as.integer(n.obs), as.double(loc.dsgn.mat), as.integer(n.loccoeff),
                    as.double(scale.dsgn.mat), as.integer(n.scalecoeff), as.double(shape.dsgn.mat),
                    as.integer(n.shapecoeff), as.double(loc.param), as.double(scale.param),
                    as.double(shape.param), as.double(cov11), as.double(cov12),
                    as.double(cov22), fit.marge, as.double(weights),
                    hess = double(n.obs * n.param * n.pairs),
                    grad = double(n.obs * n.param), PACKAGE = "SpatialExtremes")
    
    else
      std.err <- .C("wsmithgrad3d", as.double(data), as.double(distVec), as.integer(n.site),
                    as.integer(n.obs), as.double(loc.dsgn.mat), as.integer(n.loccoeff),
                    as.double(scale.dsgn.mat), as.integer(n.scalecoeff), as.double(shape.dsgn.mat),
                    as.integer(n.shapecoeff), as.double(loc.param), as.double(scale.param),
                    as.double(shape.param), as.double(cov11), as.double(cov12), as.double(cov13),
                    as.double(cov22), as.double(cov23), as.double(cov33), fit.marge,
                    as.double(weights), hess = double(n.obs * n.param * n.pairs),
                    grad = double(n.obs * n.param), PACKAGE = "SpatialExtremes")
  }
  
  grad <- matrix(std.err$grad, nrow = n.obs, ncol = n.param)
  hess <- matrix(std.err$hess, nrow = n.obs * n.pairs, ncol = n.param)

  if (iso){
    if (dist.dim == 2){
      grad[,1] <- rowSums(grad[,c(1,3)])
      grad <- grad[,-(2:3), drop = FALSE]

      hess[,1] <- rowSums(hess[,c(1,3)])
      hess <- hess[,-(2:3), drop = FALSE]
    }

    if (dist.dim == 3){
      grad[,1] <- rowSums(grad[,c(1,4,6)])
      grad <- grad[,-(2:6), drop = FALSE]

      hess[,1] <- rowSums(grad[,c(1,4,6)])
      hess <- hess[,-(2:6), drop = FALSE]
    }
  }

  n.fixed <- length(fixed.param)
  if (n.fixed > 0){
    idx <- which(param.names %in% fixed.param)
    grad <- grad[, -idx, drop = FALSE]
    hess <- hess[, -idx, drop = FALSE]
  }

  if (any(is.na(grad)))
    return(list(var.score = NA, hessian = NA, gradient = NA))

  if (std.err.type == "score"){
    var.score <- var(grad) * n.obs
    hessian <- var(hess) * n.obs * n.pairs
  }
  
  if (std.err.type == "grad"){
    var.score <- matrix(0, ncol(grad), ncol(grad))
    for (i in 1:n.obs)
      var.score <- var.score + grad[i,] %*% t(grad[i,])

    hessian <- matrix(0, ncol(hess), ncol(hess))
    for (i in 1:(n.obs * n.pairs))
      hessian <- hessian + hess[i,] %*% t(hess[i,])
  }
  
  gradient <- as.double(colSums(grad))
  
  return(list(var.score = var.score, hessian = hessian, gradient = gradient))  
}

.schlatherstderr <- function(par, data, dist, cov.mod, loc.dsgn.mat,
                             scale.dsgn.mat, shape.dsgn.mat, fit.marge,
                             std.err.type = "score", fixed.param, param.names,
                             weights){

  ##data is a matrix with each column corresponds to one location
  ##distVec is the a matrix giving the "distance vector" for each pair
  ##(1 row = 1 station)
  n.site <- ncol(data)
  n.obs <- nrow(data)
  n.pairs <- n.site * (n.site - 1) / 2
  n.param <- length(param.names)

  sill <- par["sill"]
  range <- par["range"]
  smooth <- par["smooth"]

  if (cov.mod == 5)
    ##i.e. Generalized Cauchy
    smooth2 <- par["smooth2"]

  else
    ##it won't be used anyway...
    smooth2 <- 0
  
  if (fit.marge){

    n.loccoeff <- ncol(loc.dsgn.mat)
    n.scalecoeff <- ncol(scale.dsgn.mat)
    n.shapecoeff <- ncol(shape.dsgn.mat)

    loc.idx <- which(substr(names(par), 1, 3) == "loc")
    scale.idx <- which(substr(names(par), 1, 6) == "scaleC")
    shape.idx <- which(substr(names(par), 1, 5) == "shape")
    
    loc.param <- par[loc.idx]
    scale.param <- par[scale.idx]
    shape.param <- par[shape.idx]
  }

  else {
    n.loccoeff <- 1
    n.scalecoeff <- 1
    n.shapecoeff <- 1
    
    loc.param <- 1
    scale.param <- 1
    shape.param <- 1
  }

  if (is.null(weights))
    std.err <- .C("schlatherstderr", as.integer(cov.mod), as.double(data),
                  as.double(dist), as.integer(n.site),
                  as.integer(n.obs), as.double(loc.dsgn.mat),
                  as.integer(n.loccoeff), as.double(scale.dsgn.mat),
                  as.integer(n.scalecoeff), as.double(shape.dsgn.mat),
                  as.integer(n.shapecoeff), as.double(loc.param),
                  as.double(scale.param), as.double(shape.param),
                  as.double(sill), as.double(range), as.double(smooth),
                  as.double(smooth2), fit.marge, hess = double(n.obs * n.param * n.pairs),
                  grad = double(n.obs * n.param), PACKAGE = "SpatialExtremes")

  else
    std.err <- .C("wschlatherstderr", as.integer(cov.mod), as.double(data),
                  as.double(dist), as.integer(n.site),
                  as.integer(n.obs), as.double(loc.dsgn.mat),
                  as.integer(n.loccoeff), as.double(scale.dsgn.mat),
                  as.integer(n.scalecoeff), as.double(shape.dsgn.mat),
                  as.integer(n.shapecoeff), as.double(loc.param),
                  as.double(scale.param), as.double(shape.param),
                  as.double(sill), as.double(range), as.double(smooth),
                  as.double(smooth2), fit.marge, as.double(weights),
                  hess = double(n.obs * n.param * n.pairs),
                  grad = double(n.obs * n.param), PACKAGE = "SpatialExtremes")

  grad <- matrix(std.err$grad, nrow = n.obs, ncol = n.param)
  hess <- matrix(std.err$hess, nrow = n.obs * n.pairs, ncol = n.param)

  n.fixed <- length(fixed.param)
  if (n.fixed > 0){
    idx <- which(param.names %in% fixed.param)
    grad <- grad[, -idx, drop = FALSE]
    hess <- hess[, -idx, drop = FALSE]
  }

  if (any(is.na(grad)))
    return(list(var.score = NA, hessian = NA, gradient = NA))

  if (std.err.type == "score"){
    var.score <- var(grad) * n.obs
    hessian <- var(hess) * n.obs * n.pairs
  }
  
  if (std.err.type == "grad"){
    var.score <- matrix(0, ncol(grad), ncol(grad))
    for (i in 1:n.obs)
      var.score <- var.score + grad[i,] %*% t(grad[i,])

    hessian <- matrix(0, ncol(hess), ncol(hess))
    for (i in 1:(n.obs * n.pairs))
      hessian <- hessian + hess[i,] %*% t(hess[i,])
  }

  gradient <- as.double(colSums(grad))

  return(list(var.score = var.score, hessian = hessian, gradient = gradient))
}

.schlatherindstderr <- function(par, data, dist, cov.mod, loc.dsgn.mat,
                                scale.dsgn.mat, shape.dsgn.mat, fit.marge,
                                std.err.type = "score", fixed.param, param.names,
                                weights){

  ##data is a matrix with each column corresponds to one location
  ##distVec is the a matrix giving the "distance vector" for each pair
  ##(1 row = 1 station)
  n.site <- ncol(data)
  n.obs <- nrow(data)
  n.pairs <- n.site * (n.site - 1) / 2
  n.param <- length(param.names)

  alpha <- par["alpha"]
  sill <- par["sill"]
  range <- par["range"]
  smooth <- par["smooth"]

  if (cov.mod == 5)
    ##i.e. Generalized cauchy
    smooth2 <- par["smooth2"]

  else
    ##It won't be used anyway
    smooth2 <- 0
  
  if (fit.marge){
    n.loccoeff <- ncol(loc.dsgn.mat)
    n.scalecoeff <- ncol(scale.dsgn.mat)
    n.shapecoeff <- ncol(shape.dsgn.mat)

    loc.idx <- which(substr(names(par), 1, 3) == "loc")
    scale.idx <- which(substr(names(par), 1, 6) == "scaleC")
    shape.idx <- which(substr(names(par), 1, 5) == "shape")
    
    loc.param <- par[loc.idx]
    scale.param <- par[scale.idx]
    shape.param <- par[shape.idx]
  }

  else {
    n.loccoeff <- 1
    n.scalecoeff <- 1
    n.shapecoeff <- 1
    
    loc.param <- 1
    scale.param <- 1
    shape.param <- 1
  }

  if (is.null(weights))
    std.err <- .C("schlatherindstderr", as.integer(cov.mod), as.double(data),
                  as.double(dist), as.integer(n.site),
                  as.integer(n.obs), as.double(loc.dsgn.mat),
                  as.integer(n.loccoeff), as.double(scale.dsgn.mat),
                  as.integer(n.scalecoeff), as.double(shape.dsgn.mat),
                  as.integer(n.shapecoeff), as.double(loc.param),
                  as.double(scale.param), as.double(shape.param),
                  as.double(alpha), as.double(sill), as.double(range),
                  as.double(smooth), as.double(smooth2), fit.marge,
                  hess = double(n.obs * n.param * n.pairs),
                  grad = double(n.obs * n.param),
                  PACKAGE = "SpatialExtremes")

  else
    std.err <- .C("wschlatherindstderr", as.integer(cov.mod), as.double(data),
                  as.double(dist), as.integer(n.site), as.integer(n.obs),
                  as.double(loc.dsgn.mat), as.integer(n.loccoeff),
                  as.double(scale.dsgn.mat), as.integer(n.scalecoeff),
                  as.double(shape.dsgn.mat), as.integer(n.shapecoeff),
                  as.double(loc.param), as.double(scale.param),
                  as.double(shape.param), as.double(alpha), as.double(sill),
                  as.double(range), as.double(smooth), as.double(smooth2), fit.marge,
                  as.double(weights), hess = double(n.obs * n.param * n.pairs),
                  grad = double(n.obs * n.param), PACKAGE = "SpatialExtremes")

  grad <- matrix(std.err$grad, nrow = n.obs, ncol = n.param)
  hess <- matrix(std.err$hess, nrow = n.obs * n.pairs, ncol = n.param)
  
  n.fixed <- length(fixed.param)
  if (n.fixed > 0){
    idx <- which(param.names %in% fixed.param)
    grad <- grad[, -idx, drop = FALSE]
    hess <- hess[, -idx, drop = FALSE]
  }

  if (any(is.na(grad)))
    return(list(var.score = NA, hessian = NA, gradient = NA))

  if (std.err.type == "score"){
    hessian <- var(hess) * n.obs * n.pairs
    var.score <- var(grad) * n.obs
  }
  
  if (std.err.type == "grad"){
    var.score <- matrix(0, ncol(grad), ncol(grad))
    for (i in 1:n.obs)
      var.score <- var.score + grad[i,] %*% t(grad[i,])

    hessian <- matrix(0, ncol(hess), ncol(hess))
    for (i in 1:(n.obs * n.pairs))
      hessian <- hessian + hess[i,] %*% t(hess[i,])
  }

  gradient <- as.double(colSums(grad))

  return(list(var.score = var.score, hessian = hessian, gradient))
}

.geomgaussstderr <- function(par, data, dist, cov.mod, loc.dsgn.mat,
                             scale.dsgn.mat, shape.dsgn.mat, fit.marge,
                             std.err.type = "score", fixed.param, param.names,
                             weights){

  ##data is a matrix with each column corresponds to one location
  ##distVec is the a matrix giving the "distance vector" for each pair
  ##(1 row = 1 station)
  n.site <- ncol(data)
  n.obs <- nrow(data)
  n.pairs <- n.site * (n.site - 1) / 2
  n.param <- length(param.names)

  sigma2 <- par["sigma2"]
  sill <- par["sill"]
  range <- par["range"]
  smooth <- par["smooth"]

  if (cov.mod == 5)
    ##i.e. Generalized cauchy
    smooth2 <- par["smooth2"]

  else
    ##it won't be used anyway
    smooth2 <- 0
  
  if (fit.marge){

    n.loccoeff <- ncol(loc.dsgn.mat)
    n.scalecoeff <- ncol(scale.dsgn.mat)
    n.shapecoeff <- ncol(shape.dsgn.mat)

    loc.idx <- which(substr(names(par), 1, 3) == "loc")
    scale.idx <- which(substr(names(par), 1, 6) == "scaleC")
    shape.idx <- which(substr(names(par), 1, 5) == "shape")

    loc.param <- par[loc.idx]
    scale.param <- par[scale.idx]
    shape.param <- par[shape.idx]
  }

  else {
    n.loccoeff <- 1
    n.scalecoeff <- 1
    n.shapecoeff <- 1
    
    loc.param <- 1
    scale.param <- 1
    shape.param <- 1
  }

  if (is.null(weights))
    std.err <- .C("geomgaussstderr", as.integer(cov.mod), as.double(data),
                  as.double(dist), as.integer(n.site),
                  as.integer(n.obs), as.double(loc.dsgn.mat),
                  as.integer(n.loccoeff), as.double(scale.dsgn.mat),
                  as.integer(n.scalecoeff), as.double(shape.dsgn.mat),
                  as.integer(n.shapecoeff), as.double(loc.param),
                  as.double(scale.param), as.double(shape.param),
                  as.double(sigma2), as.double(sill), as.double(range),
                  as.double(smooth), as.double(smooth2), fit.marge,
                  hess = double(n.obs * n.param * n.pairs),
                    grad = double(n.obs * n.param),
                  PACKAGE = "SpatialExtremes")
  
  else
    std.err <- .C("wgeomgaussstderr", as.integer(cov.mod), as.double(data),
                  as.double(dist), as.integer(n.site),
                  as.integer(n.obs), as.double(loc.dsgn.mat),
                  as.integer(n.loccoeff), as.double(scale.dsgn.mat),
                  as.integer(n.scalecoeff), as.double(shape.dsgn.mat),
                  as.integer(n.shapecoeff), as.double(loc.param),
                  as.double(scale.param), as.double(shape.param),
                  as.double(sigma2), as.double(sill), as.double(range),
                  as.double(smooth), as.double(smooth2), fit.marge,
                  as.double(weights), hess = double(n.obs * n.param * n.pairs),
                  grad = double(n.obs * n.param),
                  PACKAGE = "SpatialExtremes")

  grad <- matrix(std.err$grad, nrow = n.obs, ncol = n.param)
  hess <- matrix(std.err$hess, nrow = n.obs * n.pairs, ncol = n.param)

  n.fixed <- length(fixed.param)
  if (n.fixed > 0){
    idx <- which(param.names %in% fixed.param)    
    grad <- grad[, -idx, drop = FALSE]
    hess <- hess[, -idx, drop = FALSE]
  }

  if (any(is.na(grad)))
    return(list(var.score = NA, hessian = NA, gradient = NA))

  if (std.err.type == "score"){
    hessian <- var(hess) * n.obs * n.pairs
    var.score <- var(grad) * n.obs
  }
  
  if (std.err.type == "grad"){
    var.score <- matrix(0, ncol(grad), ncol(grad))
    for (i in 1:n.obs)
      var.score <- var.score + grad[i,] %*% t(grad[i,])

    hessian <- matrix(0, ncol(hess), ncol(hess))
    for (i in 1:(n.obs * n.pairs))
      hessian <- hessian + hess[i,] %*% t(hess[i,])
  }
  
  gradient <- as.double(colSums(grad))

  return(list(var.score = var.score, hessian = hessian, gradient))
}


.brownresnickstderr <- function(par, data, dist, loc.dsgn.mat, scale.dsgn.mat,
                                shape.dsgn.mat, fit.marge, std.err.type = "score",
                                fixed.param, param.names, weights){

  ##data is a matrix with each column corresponds to one location
  ##distVec is the a matrix giving the "distance vector" for each pair
  ##(1 row = 1 station)
  n.site <- ncol(data)
  n.obs <- nrow(data)
  n.pairs <- n.site * (n.site - 1) / 2
  n.param <- length(param.names)

  range <- par["range"]
  smooth <- par["smooth"]

  if (fit.marge){

    n.loccoeff <- ncol(loc.dsgn.mat)
    n.scalecoeff <- ncol(scale.dsgn.mat)
    n.shapecoeff <- ncol(shape.dsgn.mat)

    loc.idx <- which(substr(names(par), 1, 3) == "loc")
    scale.idx <- which(substr(names(par), 1, 6) == "scaleC")
    shape.idx <- which(substr(names(par), 1, 5) == "shape")

    loc.param <- par[loc.idx]
    scale.param <- par[scale.idx]
    shape.param <- par[shape.idx]
  }

  else {
    n.loccoeff <- 1
    n.scalecoeff <- 1
    n.shapecoeff <- 1
    
    loc.param <- 1
    scale.param <- 1
    shape.param <- 1
  }

  if (is.null(weights))
    std.err <- .C("brownresnickstderr", as.double(data), as.double(dist), as.integer(n.site),
                  as.integer(n.obs), as.double(loc.dsgn.mat), as.integer(n.loccoeff),
                  as.double(scale.dsgn.mat), as.integer(n.scalecoeff), as.double(shape.dsgn.mat),
                  as.integer(n.shapecoeff), as.double(loc.param), as.double(scale.param),
                  as.double(shape.param), as.double(range), as.double(smooth), fit.marge,
                  hess = double(n.obs * n.param * n.pairs), grad = double(n.obs * n.param),
                  PACKAGE = "SpatialExtremes")

  else
    std.err <- .C("wbrownresnickstderr", as.double(data), as.double(dist), as.integer(n.site),
                  as.integer(n.obs), as.double(loc.dsgn.mat), as.integer(n.loccoeff),
                  as.double(scale.dsgn.mat), as.integer(n.scalecoeff), as.double(shape.dsgn.mat),
                  as.integer(n.shapecoeff), as.double(loc.param), as.double(scale.param),
                  as.double(shape.param), as.double(range), as.double(smooth), fit.marge,
                  as.double(weights), hess = double(n.obs * n.param * n.pairs),
                  grad = double(n.obs * n.param), PACKAGE = "SpatialExtremes")

  grad <- matrix(std.err$grad, nrow = n.obs, ncol = n.param)
  hess <- matrix(std.err$hess, nrow = n.obs * n.pairs, ncol = n.param)

  n.fixed <- length(fixed.param)
  if (n.fixed > 0){
    idx <- which(param.names %in% fixed.param)
    grad <- grad[, -idx, drop = FALSE]
    hess <- hess[, -idx, drop = FALSE]
  }

  if (any(is.na(grad)))
    return(list(var.score = NA, hessian = NA, gradient = NA))

  if (std.err.type == "score"){
    hessian <- var(hess) * n.obs * n.pairs
    var.score <- var(grad) * n.obs
  }
    
  if (std.err.type == "grad"){
    var.score <- matrix(0, ncol(grad), ncol(grad))
    for (i in 1:n.obs)
      var.score <- var.score + grad[i,] %*% t(grad[i,])

    hessian <- matrix(0, ncol(hess), ncol(hess))
    for (i in 1:(n.obs * n.pairs))
      hessian <- hessian + hess[i,] %*% t(hess[i,])
  }
  
  gradient <- as.double(colSums(grad))

  return(list(var.score = var.score, hessian = hessian, gradient))
}

.spatgevstderr <- function(par, data, loc.dsgn.mat, scale.dsgn.mat,
                           shape.dsgn.mat, std.err.type = "score",
                           fixed.param, param.names){

  ##data is a matrix with each column corresponds to one location
  n.site <- ncol(data)
  n.obs <- nrow(data)
  n.param <- length(param.names)
  
  n.loccoeff <- ncol(loc.dsgn.mat)
  n.scalecoeff <- ncol(scale.dsgn.mat)
  n.shapecoeff <- ncol(shape.dsgn.mat)

  loc.idx <- which(substr(names(par), 1, 3) == "loc")
  scale.idx <- which(substr(names(par), 1, 6) == "scaleC")
  shape.idx <- which(substr(names(par), 1, 5) == "shape")

  loc.param <- par[loc.idx]
  scale.param <- par[scale.idx]
  shape.param <- par[shape.idx]
  
  std.err <- .C("spatgevstderr", as.double(data), as.integer(n.site),
                as.integer(n.obs), as.double(loc.dsgn.mat),
                as.integer(n.loccoeff), as.double(scale.dsgn.mat),
                as.integer(n.scalecoeff), as.double(shape.dsgn.mat),
                as.integer(n.shapecoeff), as.double(loc.param),
                as.double(scale.param), as.double(shape.param),
                hess = double(n.obs * n.param * n.site),
                grad = double(n.obs * n.param),
                PACKAGE = "SpatialExtremes")

  grad <- matrix(std.err$grad, nrow = n.obs, ncol = n.param)
  hess <- matrix(std.err$hess, nrow = n.obs * n.site, ncol = n.param)
  
  n.fixed <- length(fixed.param)
  if (n.fixed > 0){
    idx <- which(param.names %in% fixed.param)
    hess <- hess[, -idx, drop = FALSE]
    grad <- grad[, -idx, drop = FALSE]
  }

  if (any(is.na(grad)))
    return(list(gradient = NA, var.score = NA, hessian = NA))

  if (std.err.type == "score"){
    hessian <- var(hess) * n.obs * n.site
    var.score <- var(grad) * n.obs
  }
  
  if (std.err.type == "grad"){
    var.score <- matrix(0, ncol(grad), ncol(grad))
    for (i in 1:n.obs)
      var.score <- var.score + grad[i,] %*% t(grad[i,])

    hessian <- matrix(0, ncol(hess), ncol(hess))
    for (i in 1:(n.obs * n.site))
      hessian <- hessian + hess[i,] %*% t(hess[i,])
  }

  gradient <- as.double(colSums(grad))

  return(list(var.score = var.score, hessian = hessian, gradient = gradient))
}
