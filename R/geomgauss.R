##This file contains all the function to fit the max-stable
##characterisation of Schlather using the geometric gaussian
##representation

##This functions fits the model without any spatial structure for the
##GEV parameters. Thus, each GEV parameters are estimated at each
##location. However, if fit.marge = FALSE, observation are supposed to
##be unit Frechet and only the covariance function parameters are
##estimated.
geomgaussfull <- function(data, coord, start, cov.mod = "whitmat", ...,
                          fit.marge = FALSE, warn = TRUE, method = "BFGS",
                          control = list(), std.err.type = "none", corr = FALSE){
  ##data is a matrix with each column corresponds to one location
  ##locations is a matrix giving the coordinates (1 row = 1 station)
  n.site <- ncol(data)
  n.obs <- nrow(data)
  dist.dim <- ncol(coord)
  n.pairs <- n.site * (n.site - 1) / 2

  dist <- distance(coord)
  
  if (!(cov.mod %in% c("whitmat","cauchy","powexp", "bessel")))
    stop("''cov.mod'' must be one of 'whitmat', 'cauchy', 'powexp', 'bessel'")

  if (cov.mod == "whitmat")
    cov.mod.num <- 1
  if (cov.mod == "cauchy")
    cov.mod.num <- 2
  if (cov.mod == "powexp")
    cov.mod.num <- 3
  if (cov.mod == "bessel")
    cov.mod.num <- 4
  
  ##First create a "void" function
  nplk <- function(x) x

  ##And define the "body" of the function as the number of parameters
  ##to estimate depends on n.site
  if (fit.marge){
    loc.names <- paste("loc", 1:n.site, sep="")
    scale.names <- paste("scale", 1:n.site, sep="")
    shape.names <- paste("shape", 1:n.site, sep="")
    
    param <- c("sigma2", "sill", "range", "smooth", loc.names, scale.names, shape.names)

    body(nplk) <- parse(text = paste("-.C('geomgaussfull', as.integer(cov.mod.num), as.double(data), as.double(dist), as.integer(n.site), as.integer(n.obs), as.integer(dist.dim),",
                            paste("as.double(c(", paste(loc.names, collapse = ","), ")), "),
                            paste("as.double(c(", paste(scale.names, collapse = ","), ")), "),
                            paste("as.double(c(", paste(shape.names, collapse = ","), ")), "),
                            "as.double(sigma2), as.double(sill), as.double(range), as.double(smooth), fit.marge, dns = double(1), PACKAGE = 'SpatialExtremes')$dns"))
  }

  else{
    body(nplk) <- parse(text = paste("-.C('geomgaussfull', as.integer(cov.mod.num), as.double(data), as.double(dist), as.integer(n.site), as.integer(n.obs), as.integer(dist.dim),",
                            paste("as.double(rep(1,", n.site, ")), "),
                            paste("as.double(rep(1,", n.site, ")), "),
                            paste("as.double(rep(1,", n.site, ")), "),
                            "as.double(sigma2), as.double(sill), as.double(range), as.double(smooth), fit.marge, dns = double(1), PACKAGE = 'SpatialExtremes')$dns"))
    param <- c("sigma2", "sill", "range", "smooth")
  }

  fixed.param <- list(...)[names(list(...)) %in% param]

  if ((cov.mod == "whitmat") && !("smooth" %in% names(fixed.param)) && (std.err.type != "none")){
    if (warn)
      warning("The Whittle-Matern covariance function is not differentiable w.r.t. the ''smooth'' parameter
Standard errors are not available unless you fix it.")
    
    std.err.type <- "none"
  }

  if ((cov.mod == "bessel") && !("smooth" %in% names(fixed.param)) && (std.err.type != "none")){
    if (warn)
      warning("The Bessel covariance function is not differentiable w.r.t. the ''smooth'' parameter
Standard errors are not available unless you fix it.")
    
    std.err.type <- "none"
  }
  
  if (std.err.type == "none")
    hessian <- FALSE

  else
    hessian <- TRUE

  ##Define the formal arguments of the function
  form.nplk <- NULL
  for (i in 1:length(param))
    form.nplk <- c(form.nplk, alist(a=))

  names(form.nplk) <- param
  formals(nplk) <- form.nplk
  
  if (missing(start)) {

    start <- list()
    if (fit.marge){
      locs <- scales <- rep(NA, n.site)
      shapes <- rep(0, n.site)
      
      for (i in 1:n.site){
        marg.param <- gevmle(data[,i])
        locs[i] <- marg.param["loc"]
        scales[i] <- marg.param["scale"]
        shapes[i] <- marg.param["shape"]
      }
      
      start <- as.list(unlist(list(loc = locs, scale = scales, shape = shapes)))
    }

    if (length(fixed.param) > 0){
      if (any(names(fixed.param) == "sigma2")){
        idx <- which(names(fixed.param) == "sigma2")
        args <- c(list(data = data, coord = coord, cov.mod = cov.mod, marge = "emp"),
                  fixed.param[-idx])
      }

      else
        args <- c(list(data = data, coord = coord, cov.mod = paste("g", cov.mod, sep=""),
                       marge = "emp"), fixed.param)
      
      cov.start <- do.call("fitcovariance", args)$param
    }

    else
      cov.start <- fitcovariance(data, coord, paste("g", cov.mod, sep=""),
                                 marge = "emp")$param

    start <- c(as.list(cov.start), start)
    start <- start[!(param %in% names(list(...)))]
  }
  

  if (!is.list(start)) 
    stop("'start' must be a named list")

  if (!length(start)) 
    stop("there are no parameters left to maximize over")

  nm <- names(start)
  l <- length(nm)
  f <- formals(nplk)
  names(f) <- param
  m <- match(nm, param)
  
  if(any(is.na(m))) 
    stop("'start' specifies unknown arguments")

  formals(nplk) <- c(f[m], f[-m])
  nllh <- function(p, ...) nplk(p, ...)

  if(l > 1)
    body(nllh) <- parse(text = paste("nplk(", paste("p[",1:l,
                            "]", collapse = ", "), ", ...)"))
  
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  
  start.arg <- c(list(p = unlist(start)), fixed.param)

  init.lik <- do.call("nllh", start.arg)
  if (warn && (init.lik >= 1.0e15)) 
    warning("negative log-likelihood is infinite at starting values")

  if (method == "nlminb"){
    start <- as.numeric(start)
    opt <- nlminb(start, nllh, ..., control = control)
    opt$counts <- opt$evaluations
    opt$value <- opt$objective
    names(opt$par) <- nm
    
    if ((opt$convergence != 0) || (opt$value >= 1.0e15)) {
      if (warn)
        warning("optimization may not have succeeded")
    }

    if (opt$convergence == 0)
      opt$convergence <- "successful"
  }
  
  if (method == "nlm"){
    start <- as.numeric(start)
    opt <- nlm(nllh, start, hessian = hessian, ...)
    opt$counts <- opt$iterations
    names(opt$counts) <- "function"
    opt$value <- opt$minimum
    opt$par <- opt$estimate
    names(opt$par) <- nm

    if (opt$code <= 2)
      opt$convergence <- "sucessful"
    
    if (opt$code == 3)
      opt$convergence <- "local minimum or 'steptol' is too small"

    if (opt$code == 4)
      opt$convergence <- "iteration limit reached"

    if (opt$code == 5)
      opt$convergence <- "optimization failed"
  }

  if (!(method %in% c("nlm", "nlminb"))){
    opt <- optim(start, nllh, hessian = hessian, ..., method = method,
                 control = control)
  
    if ((opt$convergence != 0) || (opt$value >= 1.0e15)) {
      
      if (warn)
        warning("optimization may not have succeeded")
      
      if (opt$convergence == 1) 
        opt$convergence <- "iteration limit reached"
    }

    else opt$convergence <- "successful"
  }

  if (opt$value == init.lik){
    if (warn)
      warning("optimization stayed at the starting values.")
    
    opt$convergence <- "Stayed at start. val."
  }

  param.names <- param
  param <- c(opt$par, unlist(fixed.param))
  param <- param[param.names]

  if ((cov.mod == "whitmat") && !("smooth" %in% names(fixed.param)) && (std.err.type != "none")){
    if (warn)
      warning("The Whittle-Matern covariance function is not differentiable w.r.t. the ''smooth'' parameter
Standard errors are not available unless you fix it.")
    
    std.err.type <- "none"
  }

  if ((cov.mod == "bessel") && !("smooth" %in% names(fixed.param)) && (std.err.type != "none")){
    if (warn)
      warning("The Bessel covariance function is not differentiable w.r.t. the ''smooth'' parameter
Standard errors are not available unless you fix it.")
    
    std.err.type <- "none"
  }
  
  if (std.err.type != "none"){
    
    var.cov <- try(solve(opt$hessian), silent = TRUE)
    if(!is.matrix(var.cov)){
      if (warn)
        warning("observed information matrix is singular; passing std.err.type to ''none''")
      
      std.err.type <- "none"
      return
    }

    else{
      ihessian <- var.cov
      jacobian <- .geomgaussgrad(param, data, dist, cov.mod.num, as.double(0),
                                 as.double(0), as.double(0), fit.marge = fit.marge,
                                 std.err.type = std.err.type, fixed.param = names(fixed.param),
                                 param.names = param.names)

      if(any(is.na(jacobian))){
        if (warn)
          warning("observed information matrix is singular; passing std.err.type to ''none''")
        
        std.err.type <- "none"
      }
    }

    if (std.err.type != "none"){      
      var.cov <- var.cov %*% jacobian %*% var.cov
      std.err <- diag(var.cov)

      std.idx <- which(std.err <= 0)
      if(length(std.idx) > 0){
        if (warn)
          warning("Some (observed) standard errors are negative;\n passing them to NA")
        
        std.err[std.idx] <- NA
      }
      
      std.err <- sqrt(std.err)
      
      if(corr) {
        .mat <- diag(1/std.err, nrow = length(std.err))
        corr.mat <- structure(.mat %*% var.cov %*% .mat, dimnames = list(nm,nm))
        diag(corr.mat) <- rep(1, length(std.err))
      }
      
      else
        corr.mat <- NULL
      
      colnames(var.cov) <- rownames(var.cov) <- colnames(ihessian) <- 
        rownames(ihessian) <- names(std.err) <- nm
    }
  }

  if (std.err.type == "none"){
    std.err <- std.err.type <- corr.mat <- NULL
    var.cov <- ihessian <- jacobian <- NULL
  }

  cov.fun <-  covariance(sill = param["sill"], range = param["range"],
                         smooth = param["smooth"], cov.mod = cov.mod, plot = FALSE)
  
  ext.coeff <- function(h)
    2 * pnorm(sqrt(param["sigma2"] * (1 - cov.fun(h)) / 2))

  fitted <- list(fitted.values = opt$par, std.err = std.err, std.err.type = std.err.type,
                 var.cov = var.cov, param = param, cov.fun = cov.fun, fixed = unlist(fixed.param),
                 deviance = 2*opt$value, corr = corr.mat, convergence = opt$convergence,
                 counts = opt$counts, message = opt$message, est = "MPLE", data = data,
                 logLik = -opt$value, opt.value = opt$value, model = "Geometric",
                 cov.mod = cov.mod, fit.marge = fit.marge, ext.coeff = ext.coeff,
                 hessian = opt$hessian, lik.fun = nllh, coord = coord, ihessian = ihessian,
                 jacobian = jacobian, marg.cov = NULL, nllh = nllh)
  
  class(fitted) <- c(fitted$model, "maxstab")
  return(fitted)
}

##This functions fits the model from a generic R formula
##i.e. classical linear models as well as smoothing splines with
##radial basis functions may be used to model spatially the GEV
##parameters
geomgaussform <- function(data, coord, cov.mod, loc.form, scale.form, shape.form,
                          start, fit.marge = TRUE, marg.cov = NULL, ...,
                          warn = TRUE, method = "BFGS", control = list(),
                          std.err.type = "none", corr = FALSE){
  ##data is a matrix with each column corresponds to one location
  ##coord is a matrix giving the coordinates (1 row = 1 station)
  n.site <- ncol(data)
  n.obs <- nrow(data)
  dist.dim <- ncol(coord)
  n.pair <- n.site * (n.site - 1) / 2

  dist <- distance(coord)
     
  if (!(cov.mod %in% c("whitmat","cauchy","powexp", "bessel")))
    stop("''cov.mod'' must be one of 'whitmat', 'cauchy', 'powexp', 'bessel'")

  if (cov.mod == "whitmat")
    cov.mod.num <- 1
  if (cov.mod == "cauchy")
    cov.mod.num <- 2
  if (cov.mod == "powexp")
    cov.mod.num <- 3
  if (cov.mod == "bessel")
    cov.mod.num <- 4

  ##With our notation, formula must be of the form y ~ xxxx
  loc.form <- update(loc.form, y ~ .)
  scale.form <- update(scale.form, y ~ .)
  shape.form <- update(shape.form, y ~ .)

  if (is.null(marg.cov))
    covariables <- data.frame(coord)

  else
    covariables <- data.frame(coord, marg.cov)
  
  loc.model <- modeldef(covariables, loc.form)
  scale.model <- modeldef(covariables, scale.form)
  shape.model <- modeldef(covariables, shape.form)

  loc.dsgn.mat <- loc.model$dsgn.mat
  scale.dsgn.mat <- scale.model$dsgn.mat
  shape.dsgn.mat <- shape.model$dsgn.mat

  loc.pen.mat <- loc.model$pen.mat
  scale.pen.mat <- scale.model$pen.mat
  shape.pen.mat <- shape.model$pen.mat

  loc.penalty <- loc.model$penalty.tot
  scale.penalty <- scale.model$penalty.tot
  shape.penalty <- shape.model$penalty.tot

  loc.type <- loc.model$type
  scale.type <- scale.model$type
  shape.type <- shape.model$type

  ##The total number of parameters to be estimated for each GEV
  ##parameter
  n.loccoeff <- ncol(loc.dsgn.mat)
  n.scalecoeff <- ncol(scale.dsgn.mat)
  n.shapecoeff <- ncol(shape.dsgn.mat)

  ##The number of ``purely parametric'' parameters to estimate i.e. we
  ##do not consider the weigths given to each basis function
  n.pparloc <- loc.model$n.ppar
  n.pparscale <- scale.model$n.ppar
  n.pparshape <- shape.model$n.ppar
  
  loc.names <- paste("locCoeff", 1:n.loccoeff, sep="")
  scale.names <- paste("scaleCoeff", 1:n.scalecoeff, sep="")
  shape.names <- paste("shapeCoeff", 1:n.shapecoeff, sep="")
  
  param <- c("sigma2", "sill", "range", "smooth", loc.names, scale.names, shape.names)

  ##First create a "void" function
  nplk <- function(x) x

  ##And define the "body" of the function as the number of parameters
  ##to estimate depends on n.site
   body(nplk) <- parse(text = paste("-.C('geomgaussdsgnmat', as.integer(cov.mod.num), as.double(data), as.double(dist), as.integer(n.site), as.integer(n.obs), as.integer(dist.dim), as.double(loc.dsgn.mat), as.double(loc.pen.mat), as.integer(n.loccoeff), as.integer(n.pparloc), as.double(loc.penalty), as.double(scale.dsgn.mat), as.double(scale.pen.mat), as.integer(n.scalecoeff), as.integer(n.pparscale), as.double(scale.penalty), as.double(shape.dsgn.mat), as.double(shape.pen.mat), as.integer(n.shapecoeff), as.integer(n.pparshape), as.double(shape.penalty),",
                          paste("as.double(c(", paste(loc.names, collapse = ","), ")), "),
                          paste("as.double(c(", paste(scale.names, collapse = ","), ")), "),
                          paste("as.double(c(", paste(shape.names, collapse = ","), ")), "),
                         "as.double(sigma2), as.double(sill), as.double(range), as.double(smooth), dns = double(1), PACKAGE = 'SpatialExtremes')$dns"))
  ##Define the formal arguments of the function
  form.nplk <- NULL
  for (i in 1:length(param))
    form.nplk <- c(form.nplk, alist(a=))

  names(form.nplk) <- param
  formals(nplk) <- form.nplk

  if (missing(start)) {

    start <- .start.geomgauss(data, coord, covariables, cov.mod, loc.form,
                              scale.form, shape.form, method = method, ...)
    
    start <- start[!(param %in% names(list(...)))]
  
  }

  if (!is.list(start)) 
    stop("'start' must be a named list")
  
  if (!length(start)) 
    stop("there are no parameters left to maximize over")
  
  nm <- names(start)
  l <- length(nm)
  f <- formals(nplk)
  names(f) <- param
  m <- match(nm, param)
  
  if(any(is.na(m))) 
    stop("'start' specifies unknown arguments")

  formals(nplk) <- c(f[m], f[-m])
  nllh <- function(p, ...) nplk(p, ...)

  if(l > 1)
    body(nllh) <- parse(text = paste("nplk(", paste("p[",1:l,
                            "]", collapse = ", "), ", ...)"))
  
  fixed.param <- list(...)[names(list(...)) %in% param]

  if ((cov.mod == "whitmat") && !("smooth" %in% names(fixed.param)) && (std.err.type != "none")){
    if (warn)
      warning("The Whittle-Matern covariance function is not differentiable w.r.t. the ''smooth'' parameter
Standard errors are not available unless you fix it.")
    
    std.err.type <- "none"
  }

  if ((cov.mod == "bessel") && !("smooth" %in% names(fixed.param)) && (std.err.type != "none")){
    if (warn)
      warning("The Bessel covariance function is not differentiable w.r.t. the ''smooth'' parameter
Standard errors are not available unless you fix it.")
    
    std.err.type <- "none"
  }

  if (std.err.type == "none")
    hessian <- FALSE

  else
    hessian <- TRUE

  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  
  start.arg <- c(list(p = unlist(start)), fixed.param)

  init.lik <- do.call("nllh", start.arg)
  if (warn && (init.lik >= 1.0e15)) 
    warning("negative log-likelihood is infinite at starting values")

  if (method == "nlminb"){
    start <- as.numeric(start)
    opt <- nlminb(start, nllh, ..., control = control)
    opt$counts <- opt$evaluations
    opt$value <- opt$objective
    names(opt$par) <- nm
    
    if ((opt$convergence != 0) || (opt$value >= 1.0e15)) {
      if (warn)
        warning("optimization may not have succeeded")
    }

    if (opt$convergence == 0)
      opt$convergence <- "successful"
  }
  
  if (method == "nlm"){
    start <- as.numeric(start)
    opt <- nlm(nllh, start, hessian = hessian, ...)
    opt$counts <- opt$iterations
    names(opt$counts) <- "function"
    opt$value <- opt$minimum
    opt$par <- opt$estimate
    names(opt$par) <- nm

    if (opt$code <= 2)
      opt$convergence <- "sucessful"
    
    if (opt$code == 3)
      opt$convergence <- "local minimum or 'steptol' is too small"

    if (opt$code == 4)
      opt$convergence <- "iteration limit reached"

    if (opt$code == 5)
      opt$convergence <- "optimization failed"

  }

  if (!(method %in% c("nlm", "nlminb"))){
    opt <- optim(start, nllh, hessian = hessian, ..., method = method,
                 control = control)
    
    if ((opt$convergence != 0) || (opt$value >= 1.0e15)){
      if (warn)
        warning("optimization may not have succeeded")
      
      if (opt$convergence != 0) 
        opt$convergence <- "iteration limit reached"
    }

    else opt$convergence <- "successful"
  }

  if (opt$value == init.lik){
    if (warn)
      warning("optimization stayed at the starting values.")
    
    opt$convergence <- "Stayed at start. val."
  }

  param.names <- param
  param <- c(opt$par, unlist(fixed.param))
  param <- param[param.names]

  if ((cov.mod == "whitmat") && !("smooth" %in% names(fixed.param)) && (std.err.type != "none")){
    if (warn)
      warning("The Whittle-Matern covariance function is not differentiable w.r.t. the ''smooth'' parameter
Standard errors are not available unless you fix it.")
    
    std.err.type <- "none"
  }

  if ((cov.mod == "bessel") && !("smooth" %in% names(fixed.param)) && (std.err.type != "none")){
    if (warn)
      warning("The Bessel covariance function is not differentiable w.r.t. the ''smooth'' parameter
Standard errors are not available unless you fix it.")
    
    std.err.type <- "none"
  }
  
  if (std.err.type != "none"){
    
    var.cov <- try(solve(opt$hessian), silent = TRUE)
    if(!is.matrix(var.cov)){
      if (warn)
        warning("observed information matrix is singular; passing std.err.type to ''none''")
      
      std.err.type <- "none"
      return
    }

    else{
      ihessian <- var.cov
      jacobian <- .geomgaussgrad(param, data, dist, cov.mod.num, loc.dsgn.mat,
                                 scale.dsgn.mat, shape.dsgn.mat,
                                 fit.marge = fit.marge, std.err.type = std.err.type,
                                 fixed.param = names(fixed.param), param.names =
                                 param.names)

      if(any(is.na(jacobian))){
        if (warn)
          warning("observed information matrix is singular; passing std.err.type to ''none''")
        
        std.err.type <- "none"
      }
    }

    if (std.err.type != "none"){      
      var.cov <- var.cov %*% jacobian %*% var.cov
      
      std.err <- diag(var.cov)
      
      std.idx <- which(std.err <= 0)
      if(length(std.idx) > 0){
        if (warn)
          warning("Some (observed) standard errors are negative;\n passing them to NA")
        
        std.err[std.idx] <- NA
      }
      
      
      std.err <- sqrt(std.err)
      
      if(corr) {
        .mat <- diag(1/std.err, nrow = length(std.err))
        corr.mat <- structure(.mat %*% var.cov %*% .mat, dimnames = list(nm,nm))
        diag(corr.mat) <- rep(1, length(std.err))
      }
      
      else
        corr.mat <- NULL
      
      colnames(var.cov) <- rownames(var.cov) <- colnames(ihessian) <- 
        rownames(ihessian) <- names(std.err) <- nm
    }
  }

  if (std.err.type == "none"){
    std.err <- std.err.type <- corr.mat <- NULL
    var.cov <- ihessian <- jacobian <- NULL
  }

  cov.fun <- covariance(sill = param["sill"], range = param["range"],
                        smooth = param["smooth"], cov.mod = cov.mod, plot = FALSE)
  
  ext.coeff <- function(h)
    2 * pnorm(sqrt(param["sigma2"] * (1 - cov.fun(h)) / 2))
  
  fitted <- list(fitted.values = opt$par, std.err = std.err, std.err.type = std.err.type,
                 var.cov = var.cov, fixed = unlist(fixed.param), param = param,
                 deviance = 2*opt$value, corr = corr.mat, convergence = opt$convergence,
                 counts = opt$counts, message = opt$message, data = data, est = "MPLE",
                 logLik = -opt$value, opt.value = opt$value, model = "Geometric", coord = coord,
                 fit.marge = fit.marge, ext.coeff = ext.coeff, cov.mod = cov.mod, cov.fun = cov.fun,
                 loc.form = loc.form, scale.form = scale.form, shape.form = shape.form,
                 lik.fun = nllh, loc.type = loc.type, scale.type = scale.type,
                 shape.type = shape.type, ihessian = ihessian, jacobian = jacobian,
                 marg.cov = marg.cov, nllh = nllh)
  
  class(fitted) <- c(fitted$model, "maxstab")
  return(fitted)
}
