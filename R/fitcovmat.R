fitcovmat <- function(data, coord, marge = "mle", start, ...){

  n.site <- ncol(data)
  n.pairs <- n.site * (n.site - 1) / 2
  dist.dim <- ncol(coord)

  extcoeff <- fitextcoeff(data, coord, estim = "Smith",
                          plot = FALSE, loess = FALSE,
                          marge = marge)
  weights <- extcoeff[,"std.err"]
  extcoeff <- extcoeff[,"ext.coeff"]

  ##This is required otherwise we'll give to much weight this observations or because
  ##starting values could not be computed
  idx <- which((extcoeff > 1) & (extcoeff < 2) & (weights > (1e-3 * median(weights))))
  extcoeff <- extcoeff[idx]
  weights <- weights[idx]
  n.pairs.real <- length(idx)
  dist <- distance(coord, vec = TRUE)[idx,]

  if (dist.dim == 2){
    
    param <- c("cov11", "cov12", "cov22")
    
    fun2d <- function(cov11, cov12, cov22)
      .C("fitcovmat2d", as.double(cov11), as.double(cov12),
         as.double(cov22), as.integer(n.pairs.real), as.double(dist),
         as.double(extcoeff), as.double(weights), ans = double(1),
         PACKAGE = "SpatialExtremes")$ans
  }

  if (dist.dim == 3){
    param <- c("cov11", "cov12", "cov13", "cov22", "cov23", "cov33")

    fun3d <- function(cov11, cov12, cov13, cov22, cov23, cov33)
      .C("fitcovmat3d", as.double(cov11), as.double(cov12), as.double(cov13),
         as.double(cov22), as.double(cov23), as.double(cov33),
         as.integer(n.pairs.real), as.double(dist), as.double(extcoeff),
         as.double(weights), ans = double(1), PACKAGE = "SpatialExtremes")$ans
  }

  fixed.param <- list(...)[names(list(...)) %in% param]

  if (missing(start)){
    if (dist.dim == 2){
      if (any(names(fixed.param) == "cov12"))
        start <- list(cov11 = 1 + 2 * abs(list(...)$cov12),
                      cov12 = list(...)$cov12,
                      cov22 = 1 + 2 * abs(list(...)$cov12))
      
      else{
        a <- 2 * qnorm(extcoeff / 2)
        sigma.start <- mean(rowSums(dist^2) / a)
        start <- list(cov11 = sigma.start, cov12 = 0, cov22 = sigma.start)
      }
    }
    
    if (dist.dim == 3){
      a <- 2 * qnorm(extcoeff / 2)
      sigma.start <- mean(rowSums(dist^2) / a)
      start <- list(cov11 = sigma.start, cov12 = 0, cov13 = 0, cov22 = sigma.start,
                    cov23 = 0, cov33 = sigma.start)
    }
  }

  start <- start[!(param %in% names(list(...)))]

  if (!is.list(start)) 
    stop("'start' must be a named list")
  
  if (!length(start)) 
    stop("there are no parameters left to maximize over")
  
  nm <- names(start)
  l <- length(nm)

  if (dist.dim == 2)
    f <- formals(fun2d)

  if (dist.dim == 3)
    f <- formals(fun3d)
  
  names(f) <- param
  m <- match(nm, param)
  
  if(any(is.na(m))) 
    stop("'start' specifies unknown arguments")

  if (dist.dim == 2){
    formals(fun2d) <- c(f[m], f[-m])
    obj.fun <- function(p, ...) fun2d(p, ...)
    

    if (l > 1)
      body(obj.fun) <- parse(text = paste("fun2d(", paste("p[",1:l,
                               "]", collapse = ", "), ", ...)"))
  }

  if (dist.dim == 3){
    formals(fun3d) <- c(f[m], f[-m])
    obj.fun <- function(p, ...) fun3d(p, ...)
    

    if (l > 1)
      body(obj.fun) <- parse(text = paste("fun3d(", paste("p[",1:l,
                               "]", collapse = ", "), ", ...)"))
  }
    
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")

  opt <- nlm(obj.fun, unlist(start), hessian = FALSE, ...)

  if (opt$code == 4) 
    opt$convergence <- "iteration limit reached"

  if (opt$code <= 2)
    opt$convergence <- "successful"

  if (opt$code %in% c(3,5))
    opt$convergence <- "Optimization may have failed"

  param.names <- param
  names(opt$estimate) <- names(start)
  param <- c(opt$estimate, unlist(fixed.param))
  param <- param[param.names]

  if (dist.dim == 2)
    Sigma <- matrix(c(param["cov11"], param["cov12"], param["cov12"],
                      param["cov22"]), 2, 2)

  else
    Sigma <- matrix(c(param["cov11"], param["cov12"], param["cov13"],
                      param["cov12"], param["cov22"], param["cov23"],
                      param["cov13"], param["cov23"], param["cov33"]),
                    3, 3)
  
  iSigma <- solve(Sigma)

  ext.coeff <- function(posVec)
    2 * pnorm(sqrt(posVec %*% iSigma %*% posVec) / 2)

  names(opt$iterations) <- "function"
  fitted <- list(fitted.values = opt$estimate, fixed = unlist(fixed.param),
                 param = param, convergence = opt$convergence,
                 counts = opt$iterations, message = opt$message, data = data,
                 est = "Least Square", opt.value = opt$minimum, model = "Smith",
                 coord = coord, fit.marge = FALSE, cov.mod = "Gaussian",
                 ext.coeff = ext.coeff)

  class(fitted) <- c(fitted$model, "maxstab")
  return(fitted)
  
}

fitcovariance <- function(data, coord, cov.mod, marge = "mle", start,
                          ...){

  n.site <- ncol(data)
  n.pairs <- n.site * (n.site - 1) / 2
  dist.dim <- ncol(coord)

  if (!(cov.mod %in% c("whitmat","cauchy","powexp")))
    stop("''cov.mod'' must be one of 'whitmat', 'cauchy', 'powexp'")

  if (cov.mod == "whitmat")
    cov.mod.num <- 1
  if (cov.mod == "cauchy")
    cov.mod.num <- 2
  if (cov.mod == "powexp")
    cov.mod.num <- 3
  
  extcoeff <- fitextcoeff(data, coord, estim = "Smith",
                          plot = FALSE, loess = FALSE,
                          marge = marge)
  weights <- extcoeff[,"std.err"]
  extcoeff <- extcoeff[,"ext.coeff"]

  ##This is required otherwise we'll give to much weight this observations
  idx <- which(weights > (1e-3 * median(weights)))
  weights <- weights[idx]
  extcoeff <- extcoeff[idx]  
  dist <- distance(coord)[idx]
  n.pairs.real <- length(idx)

  param <- c("sill", "range", "smooth")
    
  fun <- function(sill, range, smooth)
    .C("fitcovariance", as.integer(cov.mod.num), as.double(sill), as.double(range),
       as.double(smooth), as.integer(n.pairs.real), as.double(dist),
       as.double(extcoeff), as.double(weights), ans = double(1),
       PACKAGE = "SpatialExtremes")$ans

  fixed.param <- list(...)[names(list(...)) %in% param]

  if (missing(start))
    start <- list(sill = .5, range = 0.75 * max(dist), smooth = .5)
  
  start <- start[!(param %in% names(list(...)))]

  if (!is.list(start)) 
    stop("'start' must be a named list")
  
  if (!length(start)) 
    stop("there are no parameters left to maximize over")
  
  nm <- names(start)
  l <- length(nm)
  
  f <- formals(fun)
  names(f) <- param
  m <- match(nm, param)
  
  if(any(is.na(m))) 
    stop("'start' specifies unknown arguments")

  formals(fun) <- c(f[m], f[-m])
  obj.fun <- function(p, ...) fun(p, ...)
    
  if (l > 1)
    body(obj.fun) <- parse(text = paste("fun(", paste("p[",1:l,
                             "]", collapse = ", "), ", ...)"))
  
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")

  opt <- nlm(obj.fun, unlist(start), hessian = FALSE, ...)

  if (opt$code == 4) 
    opt$convergence <- "iteration limit reached"

  if (opt$code <= 2)
    opt$convergence <- "successful"

  if (opt$code %in% c(3,5))
    opt$convergence <- "Optimization may have failed"

  param.names <- param
  names(opt$estimate) <- names(start)
  param <- c(opt$estimate, unlist(fixed.param))
  param <- param[param.names]

  cov.fun <- covariance(sill = param["sill"], range = param["range"],
                        smooth = param["smooth"], cov.mod = cov.mod, plot = FALSE)
  
  ext.coeff <- function(h)
    1 + sqrt(1 - 1/2 * (cov.fun(h) + 1))

  names(opt$iterations) <- "function"
  fitted <- list(fitted.values = opt$estimate, fixed = unlist(fixed.param),
                 param = param, convergence = opt$convergence,
                 counts = opt$iterations, message = opt$message, data = data,
                 est = "Least Square", opt.value = opt$minimum, model = "Schlather",
                 coord = coord, fit.marge = FALSE, cov.mod = cov.mod,
                 cov.fun = cov.fun, ext.coeff = ext.coeff)

  class(fitted) <- c(fitted$model, "maxstab")
  return(fitted)
  
}
