fitcovmat <- function(data, coord, marge = "mle", iso = FALSE, ..., start){

  n.site <- ncol(data)
  n.pairs <- n.site * (n.site - 1) / 2
  dist.dim <- ncol(coord)

  extcoeff <- fitextcoeff(data, coord, estim = "Smith",
                          plot = FALSE, loess = FALSE,
                          marge = marge)
  weights <- extcoeff[,"std.err"]
  extcoeff <- extcoeff[,"ext.coeff"]
  dist <- distance(coord, vec = TRUE)

  ##Check if there are really small weights, this could happen in few
  ##cases
  weights[weights <= 1e-4] <- mean(weights)

  if (dist.dim == 2){

    if (iso){
      param <- "cov"
      
      fun2diso <- function(cov)
        .C("fitcovmat2d", as.double(cov), as.double(0.0),
           as.double(cov), as.integer(n.pairs), as.double(dist),
           as.double(extcoeff), as.double(weights), ans = double(1),
           PACKAGE = "SpatialExtremes")$ans
    }

    else{
      param <- c("cov11", "cov12", "cov22")
      
      fun2d <- function(cov11, cov12, cov22)
        .C("fitcovmat2d", as.double(cov11), as.double(cov12),
           as.double(cov22), as.integer(n.pairs), as.double(dist),
           as.double(extcoeff), as.double(weights), ans = double(1),
           PACKAGE = "SpatialExtremes")$ans
    }    
  }

  if (dist.dim == 3){

    if (iso){
      param <- "cov"
      
      fun3diso <- function(cov)
        .C("fitcovmat3d", as.double(cov), as.double(0.0), as.double(0.0),
           as.double(cov), as.double(0.0), as.double(cov),
           as.integer(n.pairs), as.double(dist), as.double(extcoeff),
           as.double(weights), ans = double(1), PACKAGE = "SpatialExtremes")$ans
    }

    else{
      param <- c("cov11", "cov12", "cov13", "cov22", "cov23", "cov33")
      
      fun3d <- function(cov11, cov12, cov13, cov22, cov23, cov33)
        .C("fitcovmat3d", as.double(cov11), as.double(cov12), as.double(cov13),
           as.double(cov22), as.double(cov23), as.double(cov33),
           as.integer(n.pairs), as.double(dist), as.double(extcoeff),
           as.double(weights), ans = double(1), PACKAGE = "SpatialExtremes")$ans
    }
  }

  fixed.param <- list(...)[names(list(...)) %in% param]

  if (missing(start)){
    if (iso){
      a <- 4 * qnorm(pmin(extcoeff, 2) / 2)^2
      sigma.start <- mean(rowSums(dist^2) / a)
      start <- list(cov = sigma.start)
    }

    else{
      if (dist.dim == 2){
        if (any(names(fixed.param) == "cov12"))
          start <- list(cov11 = 1 + 2 * abs(list(...)$cov12),
                        cov12 = list(...)$cov12,
                        cov22 = 1 + 2 * abs(list(...)$cov12))
        
        else{
          a <- 4 * qnorm(pmin(extcoeff, 2) / 2)^2
          sigma.start <- mean(rowSums(dist^2) / a)
          start <- list(cov11 = sigma.start, cov12 = 0, cov22 = sigma.start)
        }
      }
    
      if (dist.dim == 3){
        a <- 4 * qnorm(pmin(extcoeff, 2) / 2)^2
        sigma.start <- mean(rowSums(dist^2) / a)
        start <- list(cov11 = sigma.start, cov12 = 0, cov13 = 0, cov22 = sigma.start,
                      cov23 = 0, cov33 = sigma.start)
      }
    }
    start <- start[!(param %in% names(list(...)))]
  }

  if (!is.list(start)) 
    stop("'start' must be a named list")
  
  if (!length(start)) 
    stop("there are no parameters left to maximize over")
  
  nm <- names(start)
  l <- length(nm)

  if (dist.dim == 2){
    if (iso)
      f <- formals(fun2diso)

    else
      f <- formals(fun2d)
  }

  if (dist.dim == 3){
    if (iso)
      f <- formals(fun3diso)

    else
      f <- formals(fun3d)
  }
  
  names(f) <- param
  m <- match(nm, param)
  
  if(any(is.na(m))) 
    stop("'start' specifies unknown arguments")

  if (dist.dim == 2){
    if (iso){
      formals(fun2diso) <- c(f[m], f[-m])
      obj.fun <- function(p, ...) fun2diso(p, ...)
    

      if (l > 1)
        body(obj.fun) <- parse(text = paste("fun2diso(", paste("p[",1:l,
                                 "]", collapse = ", "), ", ...)"))
    }

    else{
      formals(fun2d) <- c(f[m], f[-m])
      obj.fun <- function(p, ...) fun2d(p, ...)
    

      if (l > 1)
        body(obj.fun) <- parse(text = paste("fun2d(", paste("p[",1:l,
                                 "]", collapse = ", "), ", ...)"))
    }
  }

  if (dist.dim == 3){
    if (iso){
      formals(fun3diso) <- c(f[m], f[-m])
      obj.fun <- function(p, ...) fun3diso(p, ...)
      

      if (l > 1)
        body(obj.fun) <- parse(text = paste("fun3diso(", paste("p[",1:l,
                                 "]", collapse = ", "), ", ...)"))
    }

    else{
      formals(fun3d) <- c(f[m], f[-m])
      obj.fun <- function(p, ...) fun3d(p, ...)
      

      if (l > 1)
        body(obj.fun) <- parse(text = paste("fun3d(", paste("p[",1:l,
                                 "]", collapse = ", "), ", ...)"))
    }
  }
    
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")

  opt <- optim(unlist(start), obj.fun, hessian = FALSE, ...)

  if (opt$convergence == 0)
    opt$convergence <- "successful"

  else if (opt$convergence == 1) 
    opt$convergence <- "iteration limit reached"

  else
    opt$convergence <- "Optimization may have failed"

  param.names <- param
  names(opt$par) <- names(start)
  param <- c(opt$par, unlist(fixed.param))
  param <- param[param.names]

  if (iso){
    if (dist.dim == 2){
      param <- c(param["cov"], 0, param["cov"], param[-1])
      names(param)[1:3] <- c("cov11", "cov12", "cov22")
    }

    else{
      param <- c(param["cov"], 0, 0, param["cov"], 0, param["cov"],
                 param[-1])
      names(param)[1:6] <- c("cov11", "cov12", "cov13", "cov22", "cov23",
                             "cov33")
    }
  }

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

  fitted <- list(fitted.values = opt$par, fixed = unlist(fixed.param),
                 param = param, convergence = opt$convergence,
                 counts = opt$counts, message = opt$message, data = data,
                 est = "Least Square", opt.value = opt$value, model = "Smith",
                 coord = coord, fit.marge = FALSE, cov.mod = "Gaussian",
                 ext.coeff = ext.coeff)

  class(fitted) <- c(fitted$model, "maxstab")
  return(fitted)
  
}

fitcovariance <- function(data, coord, cov.mod, marge = "mle", ..., start){

  n.site <- ncol(data)
  n.pairs <- n.site * (n.site - 1) / 2
  dist.dim <- ncol(coord)

  if (substr(cov.mod, 1, 1) %in% c("i", "g")){
        
    if (substr(cov.mod, 1, 1) == "i")
      model <- "iSchlather"

    if (substr(cov.mod, 1, 1) == "g")
      model <- "Geometric"

    cov.mod <- substr(cov.mod, 2, 8)
  }

  else
    model <- "Schlather"

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
  dist <- distance(coord)

  ##Check if there are really small weights, this could happen in few
  ##cases
  weights[weights <= 1e-4] <- mean(weights)

  param <- c("sill", "range", "smooth")

  if (model == "iSchlather")
    param <- c("alpha", param)

  if (model == "Geometric")
    param <- c("sigma2", param)

  if (model == "Schlather")
    funS <- function(sill, range, smooth)
      .C("fitcovariance", as.integer(cov.mod.num), as.double(sill), as.double(range),
         as.double(smooth), as.integer(n.pairs), as.double(dist),
         as.double(extcoeff), as.double(weights), ans = double(1),
         PACKAGE = "SpatialExtremes")$ans

  else if (model == "iSchlather")
    funI <- function(alpha, sill, range, smooth)
      .C("fiticovariance", as.integer(cov.mod.num), as.double(alpha), as.double(sill),
         as.double(range), as.double(smooth), as.integer(n.pairs), as.double(dist),
         as.double(extcoeff), as.double(weights), ans = double(1),
         PACKAGE = "SpatialExtremes")$ans

  else
    funG <- function(sigma2, sill, range, smooth)
      .C("fitgcovariance", as.integer(cov.mod.num), as.double(sigma2), as.double(sill),
         as.double(range), as.double(smooth), as.integer(n.pairs), as.double(dist),
         as.double(extcoeff), as.double(weights), ans = double(1),
         PACKAGE = "SpatialExtremes")$ans

  fixed.param <- list(...)[names(list(...)) %in% param]

  if (missing(start)){
    start <- list(sill = .9, range = 0.75 * max(dist), smooth = .5)

    if (model == "iSchlather")
      start <- c(list(alpha = 0.5), start)

    if (model == "Geometric")
      start <- c(list(sigma2 = 1), start)

    start <- start[!(param %in% names(list(...)))]
  }

  if (!is.list(start)) 
    stop("'start' must be a named list")
  
  if (!length(start)) 
    stop("there are no parameters left to maximize over")
  
  nm <- names(start)
  l <- length(nm)

  if (model == "Schlather")
    f <- formals(funS)

  else if (model == "iSchlather")
    f <- formals(funI)

  else
    f <- formals(funG)
  
  names(f) <- param
  m <- match(nm, param)
  
  if(any(is.na(m))) 
    stop("'start' specifies unknown arguments")

  if (model == "Schlather"){
    formals(funS) <- c(f[m], f[-m])
    obj.fun <- function(p, ...) funS(p, ...)
    
    if (l > 1)
      body(obj.fun) <- parse(text = paste("funS(", paste("p[",1:l,
                               "]", collapse = ", "), ", ...)"))
  }

  else if (model == "iSchlather"){
    formals(funI) <- c(f[m], f[-m])
    obj.fun <- function(p, ...) funI(p, ...)
    
    if (l > 1)
      body(obj.fun) <- parse(text = paste("funI(", paste("p[",1:l,
                               "]", collapse = ", "), ", ...)"))
  }

  else {
    formals(funG) <- c(f[m], f[-m])
    obj.fun <- function(p, ...) funG(p, ...)
    
    if (l > 1)
      body(obj.fun) <- parse(text = paste("funG(", paste("p[",1:l,
                               "]", collapse = ", "), ", ...)"))
  }
  
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")

  opt <- optim(unlist(start), obj.fun, hessian = FALSE, ...)

  if (opt$convergence == 1) 
    opt$convergence <- "iteration limit reached"

  else if (opt$convergence == 0)
    opt$convergence <- "successful"

  else
    opt$convergence <- "Optimization may have failed"

  param.names <- param
  names(opt$par) <- names(start)
  param <- c(opt$par, unlist(fixed.param))
  param <- param[param.names]

  cov.fun <- covariance(sill = param["sill"], range = param["range"],
                        smooth = param["smooth"], cov.mod = cov.mod, plot = FALSE)

  if (model == "Schlather")
    ext.coeff <- function(h)
      1 + sqrt(1 - 1/2 * (cov.fun(h) + 1))

  else if (model == "iSchlather")
    ext.coeff <- function(h)
      2 * param["alpha"] + (1 - param["alpha"]) *
        (1 + sqrt(1 - 1/2 * (cov.fun(h) + 1)))

  else
    ext.coeff <- function(h)
      2 * pnorm(sqrt(param["sigma2"] * (1 - cov.fun(h)) / 2))

  fitted <- list(fitted.values = opt$par, fixed = unlist(fixed.param),
                 param = param, convergence = opt$convergence,
                 counts = opt$counts, message = opt$message, data = data,
                 est = "Least Square", opt.value = opt$value, model = model,
                 coord = coord, fit.marge = FALSE, cov.mod = cov.mod,
                 cov.fun = cov.fun, ext.coeff = ext.coeff)

  class(fitted) <- c(fitted$model, "maxstab")
  return(fitted)
  
}
