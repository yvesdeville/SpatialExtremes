distance <- function(coord, vec = FALSE){
  ##This function computes the distance between each pair of locations

  ##coord is a matrix giving the coordinates (1 row = 1 station)
  n.site <- nrow(coord)
  dist.dim <- ncol(coord)
  n.pairs <- n.site * (n.site - 1) / 2

  
  if (vec){
    dist <- .C("distance", as.double(coord), as.integer(dist.dim),
               as.integer(n.site), vec, dist = double(dist.dim * n.pairs),
               PACKAGE = "SpatialExtremes")$dist
    dist <- matrix(dist, ncol = dist.dim, nrow = n.pairs)
  }

  else    
    dist <- .C("distance", as.double(coord), as.integer(dist.dim),
               as.integer(n.site), vec, dist = double(n.pairs),
               PACKAGE = "SpatialExtremes")$dist
  
  return(dist)
}

gev2frech <- function(x, loc, scale, shape, emp = FALSE){

  if (emp){
    probs <- ppoints(x)
    x[order(x)] <- - 1 / log(probs)
    return(x)
  }
  
  if (shape == 0)
    exp((x - loc)/scale)
  
  else
    pmax(1 + shape * (x - loc) / scale, 0)^(1/shape)
}

frech2gev <- function(x, loc, scale, shape){
  if (shape == 0)
    scale * log(pmax(x, 0)) + loc

  else
    loc + scale * (pmax(x, 0)^shape - 1) / shape
}

.qgev <- function(p, loc = 1, scale = 1, shape = 1,
                  lower.tail = TRUE){
  
    if ((min(p, na.rm = TRUE) <= 0) || (max(p, na.rm = TRUE) >=1))
      stop("'p' must contain probabilities in (0,1)")
    
    if (min(scale) < 0)
      warning("There are some invalid scale GEV parameters")
    
    if (length(p) != 1)
      stop("invalid p")
    
    if (!lower.tail)
      p <- 1 - p

    n <- length(loc)

    ans <- .C("gev", as.double(p), as.integer(n), as.double(loc),
              as.double(scale), as.double(shape), quant = double(n),
              PACKAGE = "SpatialExtremes")$quant
    
    return(ans)
}
