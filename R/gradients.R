.smithgrad <- function(par, data, distVec, loc.dsgn.mat,
                       scale.dsgn.mat, shape.dsgn.mat, fit.marge,
                       std.err.type = "score", fixed.param, param.names,
                       jacobian = TRUE, iso = TRUE){

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

  if (dist.dim == 2)
    grad <- .C("smithgrad", as.double(data), as.double(distVec), as.integer(n.site),
               as.integer(n.obs), as.double(loc.dsgn.mat), as.integer(n.loccoeff),
               as.double(scale.dsgn.mat), as.integer(n.scalecoeff), as.double(shape.dsgn.mat),
               as.integer(n.shapecoeff), as.double(loc.param), as.double(scale.param),
               as.double(shape.param), as.double(cov11), as.double(cov12),
               as.double(cov22), fit.marge, grad = double(n.obs * n.param),
               PACKAGE = "SpatialExtremes")$grad

  else
    grad <- .C("smithgrad3d", as.double(data), as.double(distVec), as.integer(n.site),
               as.integer(n.obs), as.double(loc.dsgn.mat), as.integer(n.loccoeff),
               as.double(scale.dsgn.mat), as.integer(n.scalecoeff), as.double(shape.dsgn.mat),
               as.integer(n.shapecoeff), as.double(loc.param), as.double(scale.param),
               as.double(shape.param), as.double(cov11), as.double(cov12), as.double(cov13),
               as.double(cov22), as.double(cov23), as.double(cov33), fit.marge,
               grad = double(n.obs * n.param), PACKAGE = "SpatialExtremes")$grad
  
  grad <- matrix(grad, nrow = n.obs, ncol = n.param)

  if (iso){
    if (dist.dim == 2){
      grad[,1] <- rowSums(grad[,c(1,3)])
      grad <- grad[,-(2:3)]
    }

    if (dist.dim == 3){
      grad[,1] <- rowSums(grad[,c(1,4,6)])
      grad <- grad[,-(2:6)]
    }
  }

  n.fixed <- length(fixed.param)
  if (n.fixed > 0){
    idx <- NULL
    for (i in 1:n.fixed)
      idx <- c(idx, which(param.names == fixed.param[i]))
    
    grad <- grad[,-idx]
  }

  if (any(is.na(grad)))
    return(NA)

  if (jacobian){
    if (std.err.type == "score")
      jacobian <- var(grad) * n.obs
    
    if (std.err.type == "grad"){
      gradient <- colSums(grad)
      jacobian <- 0
      for (i in 1:n.obs){
        grad.vec <- matrix(grad[i,], ncol = 1)
        jacobian <- jacobian + grad.vec %*% t(grad.vec)
      }
    }
    
    return(jacobian)
  }

  else{
    return(as.double(colSums(grad)))
  }
}

.schlathergrad <- function(par, data, dist, cov.mod, loc.dsgn.mat,
                           scale.dsgn.mat, shape.dsgn.mat, fit.marge,
                           std.err.type = "score", fixed.param, param.names){

  ##data is a matrix with each column corresponds to one location
  ##distVec is the a matrix giving the "distance vector" for each pair
  ##(1 row = 1 station)
  n.site <- ncol(data)
  n.obs <- nrow(data)
  n.pairs <- n.site * (n.site - 1) / 2
    
  if (fit.marge){

    n.loccoeff <- ncol(loc.dsgn.mat)
    n.scalecoeff <- ncol(scale.dsgn.mat)
    n.shapecoeff <- ncol(shape.dsgn.mat)

    loc.idx <- which(substr(names(par), 1, 3) == "loc")
    scale.idx <- which(substr(names(par), 1, 6) == "scaleC")
    shape.idx <- which(substr(names(par), 1, 5) == "shape")

    sill <- par["sill"]
    range <- par["range"]
    smooth <- par["smooth"]
    loc.param <- par[loc.idx]
    scale.param <- par[scale.idx]
    shape.param <- par[shape.idx]
  }

  else {
    n.loccoeff <- 1
    n.scalecoeff <- 1
    n.shapecoeff <- 1
    sill <- par["sill"]
    range <- par["range"]
    smooth <- par["smooth"]
    loc.param <- 1
    scale.param <- 1
    shape.param <- 1
  }
  
  grad <- .C("schlathergrad", as.integer(cov.mod), as.double(data),
             as.double(dist), as.integer(n.site),
             as.integer(n.obs), as.double(loc.dsgn.mat),
             as.integer(n.loccoeff), as.double(scale.dsgn.mat),
             as.integer(n.scalecoeff), as.double(shape.dsgn.mat),
             as.integer(n.shapecoeff), as.double(loc.param),
             as.double(scale.param), as.double(shape.param),
             as.double(sill), as.double(range), as.double(smooth),
             fit.marge, grad = double(n.obs * length(param.names)),
             PACKAGE = "SpatialExtremes")$grad

  grad <- matrix(grad, nrow = n.obs, ncol = length(param.names))

  n.fixed <- length(fixed.param)
  if (n.fixed > 0){
    idx <- NULL
    for (i in 1:n.fixed)
      idx <- c(idx, which(param.names == fixed.param[i]))
    
    grad <- grad[,-idx]
  }

  if (any(is.na(grad)))
    return(NA)
  
  if (std.err.type == "score")
    jacobian <- var(grad) * n.obs

  if (std.err.type == "grad"){
    jacobian <- 0
    for (i in 1:n.obs){
      grad.vec <- matrix(grad[i,], ncol = 1)
      jacobian <- jacobian + grad.vec %*% t(grad.vec)
    }
  }

  return(jacobian)
}

.schlatherindgrad <- function(par, data, dist, cov.mod, loc.dsgn.mat,
                              scale.dsgn.mat, shape.dsgn.mat, fit.marge,
                              std.err.type = "score", fixed.param, param.names){

  ##data is a matrix with each column corresponds to one location
  ##distVec is the a matrix giving the "distance vector" for each pair
  ##(1 row = 1 station)
  n.site <- ncol(data)
  n.obs <- nrow(data)
  n.pairs <- n.site * (n.site - 1) / 2
    
  if (fit.marge){

    n.loccoeff <- ncol(loc.dsgn.mat)
    n.scalecoeff <- ncol(scale.dsgn.mat)
    n.shapecoeff <- ncol(shape.dsgn.mat)

    loc.idx <- which(substr(names(par), 1, 3) == "loc")
    scale.idx <- which(substr(names(par), 1, 6) == "scaleC")
    shape.idx <- which(substr(names(par), 1, 5) == "shape")

    alpha <- par["alpha"]
    sill <- par["sill"]
    range <- par["range"]
    smooth <- par["smooth"]
    loc.param <- par[loc.idx]
    scale.param <- par[scale.idx]
    shape.param <- par[shape.idx]
  }

  else {
    n.loccoeff <- 1
    n.scalecoeff <- 1
    n.shapecoeff <- 1
    alpha <- par["alpha"]
    sill <- par["sill"]
    range <- par["range"]
    smooth <- par["smooth"]
    loc.param <- 1
    scale.param <- 1
    shape.param <- 1
  }
  
  grad <- .C("schlatherindgrad", as.integer(cov.mod), as.double(data),
             as.double(dist), as.integer(n.site),
             as.integer(n.obs), as.double(loc.dsgn.mat),
             as.integer(n.loccoeff), as.double(scale.dsgn.mat),
             as.integer(n.scalecoeff), as.double(shape.dsgn.mat),
             as.integer(n.shapecoeff), as.double(loc.param),
             as.double(scale.param), as.double(shape.param),
             as.double(alpha), as.double(sill), as.double(range),
             as.double(smooth), fit.marge, grad = double(n.obs * length(param.names)),
             PACKAGE = "SpatialExtremes")$grad

  grad <- matrix(grad, nrow = n.obs, ncol = length(param.names))

  n.fixed <- length(fixed.param)
  if (n.fixed > 0){
    idx <- NULL
    for (i in 1:n.fixed)
      idx <- c(idx, which(param.names == fixed.param[i]))
    
    grad <- grad[,-idx]
  }

  if (any(is.na(grad)))
    return(NA)
  
  if (std.err.type == "score")
    jacobian <- var(grad) * n.obs

  if (std.err.type == "grad"){
    jacobian <- 0
    for (i in 1:n.obs){
      grad.vec <- matrix(grad[i,], ncol = 1)
      jacobian <- jacobian + grad.vec %*% t(grad.vec)
    }
  }

  return(jacobian)
}
