fitspatgev <- function(data, covariables, loc.form, scale.form, shape.form,
                       ..., start, control = list(maxit = 10000),
                       method = "Nelder", std.err.type = "score", warn = TRUE){

  if (std.err.type != "none")
    hessian <- std.err.flag <- TRUE

  else
    hessian <- std.err.flag <- FALSE
  
  n.site <- ncol(data)
  n.obs <- nrow(data)

  ##With our notation, formula must be of the form y ~ xxxx
  loc.form <- update(loc.form, y ~ .)
  scale.form <- update(scale.form, y ~ .)
  shape.form <- update(shape.form, y ~ .)

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

  param <- c(loc.names, scale.names, shape.names)

  nllik <- function(x) x

  body(nllik) <- parse(text = paste("-.C('spatgevlik', as.double(data), as.double(covariables), as.integer(n.site), as.integer(n.obs), as.double(loc.dsgn.mat), as.double(loc.pen.mat), as.integer(n.loccoeff), as.integer(n.pparloc), as.double(loc.penalty), as.double(scale.dsgn.mat), as.double(scale.pen.mat), as.integer(n.scalecoeff), as.integer(n.pparscale), as.double(scale.penalty), as.double(shape.dsgn.mat), as.double(shape.pen.mat), as.integer(n.shapecoeff), as.integer(n.pparshape), as.double(shape.penalty),",
                         paste("as.double(c(", paste(loc.names, collapse = ","), ")), "),
                         paste("as.double(c(", paste(scale.names, collapse = ","), ")), "),
                         paste("as.double(c(", paste(shape.names, collapse = ","), ")), "),
                         "dns = double(1), PACKAGE = 'SpatialExtremes')$dns"))

  
  ##Define the formal arguments of the function
  form.nllik <- NULL
  for (i in 1:length(param))
    form.nllik <- c(form.nllik, alist(a=))

  names(form.nllik) <- param
  formals(nllik) <- form.nllik

  if (missing(start)){
    loc <- scale <- shape <- rep(0, n.site)

    for (i in 1:n.site){
      gev.param <- gevmle(data[,i])
      loc[i] <- gev.param["loc"]
      scale[i] <- gev.param["scale"]
      shape[i] <- gev.param["shape"]
    }

    locCoeff <- loc.model$init.fun(loc)
    scaleCoeff <- scale.model$init.fun(scale)
    shapeCoeff <- shape.model$init.fun(shape)

    locCoeff[is.na(locCoeff)] <- 0
    scaleCoeff[is.na(scaleCoeff)] <- 0
    shapeCoeff[is.na(shapeCoeff)] <- 0

    ##To be sure that the scale parameter is always positive at starting
    ##values
    scales.hat <- scale.model$dsgn.mat %*% scaleCoeff
  
    if (any(scales.hat <= 0))
      scaleCoeff[1] <- scaleCoeff[1] - 1.001 * min(scales.hat)
  
    names(locCoeff) <- loc.names
    names(scaleCoeff) <- scale.names
    names(shapeCoeff) <- shape.names
    
    start <- as.list(c(locCoeff, scaleCoeff, shapeCoeff))
    start <- start[!(param %in% names(list(...)))]

  }

  if (!length(start)) 
    stop("there are no parameters left to maximize over")
  
  nm <- names(start)
  l <- length(nm)
  f <- formals(nllik)
  names(f) <- param
  m <- match(nm, param)
  
  if(any(is.na(m))) 
    stop("'start' specifies unknown arguments")

  formals(nllik) <- c(f[m], f[-m])
  nllh <- function(p, ...) nllik(p, ...)

  if(l > 1)
    body(nllh) <- parse(text = paste("nllik(", paste("p[",1:l,
                            "]", collapse = ", "), ", ...)"))
  
  fixed.param <- list(...)[names(list(...)) %in% param]
  
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  
  start.arg <- c(list(p = unlist(start)), fixed.param)

  init.lik <- do.call("nllh", start.arg)
  if (warn && (init.lik >= 1.0e6)) 
    warning("negative log-likelihood is infinite at starting values")

  if (method == "nlminb"){
    start <- as.numeric(start)
    opt <- nlminb(start, nllh, ..., control = control)
    opt$counts <- opt$evaluations
    opt$value <- opt$objective
    names(opt$par) <- nm
  }
  
  else if (method == "nlm"){
    start <- as.numeric(start)
    opt <- nlm(nllh, start, hessian = hessian, ...)
    opt$counts <- opt$iterations
    names(opt$counts) <- "function"
    opt$value <- opt$minimum
    opt$par <- opt$estimate
    names(opt$par) <- nm

    if (opt$code <= 2){
      opt$convergence <- 0
      opt$message <- NULL
    }
    
    if (opt$code > 2){
      opt$convergence <- 1
      opt$message <- paste("nlm error code", opt$code)
    }      
  }

  else 
    opt <- optim(start, nllh, hessian = hessian, ..., method = method,
                 control = control)
    
  if ((opt$convergence != 0) || (opt$value >= 1.0e6)){
    if (warn)
      warning("optimization may not have succeeded")
  }
  
  else
    opt$convergence <- "successful"

  param.names <- param
  param <- c(opt$par, unlist(fixed.param))
  param <- param[param.names]

  if (std.err.flag){
    var.cov <- try(solve(opt$hessian), silent = TRUE)
    
    if(!is.matrix(var.cov)){
      std.err.flag <- FALSE
      warning("observed information matrix is singular; std. err. won't be computed")
    }
    
    else{
      ihessian <- var.cov
      jacobian <- .spatgevgrad(param, data, loc.dsgn.mat,scale.dsgn.mat,
                               shape.dsgn.mat, std.err.type,
                               fixed.param = names(fixed.param), param.names = param.names)
      
      if (any(is.na(jacobian))){
        if (warn)
          warning("observed information matrix is singular; std. err. won't be computed")
        
        std.err.flag <- FALSE
      }
    }
  }

  if (std.err.flag){
    var.cov <- var.cov %*% jacobian %*% var.cov
    std.err <- sqrt(diag(var.cov))

    colnames(var.cov) <- colnames(ihessian) <- rownames(var.cov) <-
      rownames(ihessian) <- colnames(jacobian) <- rownames(jacobian) <-
        names(std.err) <- nm
  }

  else
    std.err <- std.err.type <- corr.mat <- var.cov <- ihessian <-
      jacobian <- NULL
  
  ans <- list(fitted.values = opt$par, param = param, std.err = std.err, var.cov = var.cov,
              counts = opt$counts, message = opt$message, covariables = covariables,
              logLik = -opt$value, loc.form = loc.form, scale.form = scale.form,
              shape.form = shape.form, convergence = opt$convergence, nllh = nllh,
              deviance = 2 * opt$value, ihessian = ihessian, jacobian = jacobian,
              data = data, jacobian = jacobian, fixed = unlist(fixed.param),
              hessian = opt$hessian)

  class(ans) <- "spatgev"
  return(ans)
}
