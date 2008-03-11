.start.smith <- function(data, coord, loc.model, scale.model, shape.model,
                         print.start.values = TRUE, method = "Nelder",
                         ...){

  if (ncol(coord) == 2)
    param <- c("cov11", "cov12", "cov22")

  else
    param <- c("cov11", "cov12", "cov13", "cov22", "cov23", "cov33")

  fixed.param <- list(...)[names(list(...)) %in% param]
  
  if (print.start.values)
    cat("Computing appropriate starting values\n")
  
  n.site <- ncol(data)

  loc <- scale <- shape <- rep(NA, n.site)
  dataFrech <- data
  for (i in 1:n.site){
    marg.param <- gevmle(data[,i])
    loc[i] <- marg.param["loc"]
    scale[i] <- marg.param["scale"]
    shape[i] <- marg.param["shape"]
    dataFrech[,i] <- gev2frech(dataFrech[,i], loc[i], scale[i], shape[i])
  }

  if (length(fixed.param) > 0){
    args <- c(list(data = data, coord = coord, marge = "emp"),
              fixed.param)
    covs <- do.call("fitcovmat", args)$param
  }

  else
    covs <- fitcovmat(data, coord, marge = "emp")$param  

  covs <- smithfull(dataFrech, coord, start = as.list(covs),
                    fit.marge = FALSE, method = method, warn = FALSE)$param
  
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
  
  start <- as.list(covs)

  names(locCoeff) <- names(scaleCoeff) <- names(shapeCoeff) <- NULL
  
  start <- c(start, as.list(unlist(list(locCoeff = locCoeff,
                                        scaleCoeff = scaleCoeff,
                                        shapeCoeff = shapeCoeff))))
  if (print.start.values)
    cat("Starting values are defined\n")

  if (print.start.values){
    cat("Starting values are:\n")
    print(unlist(start))
  }

  return(start)
}


.start.schlather <- function(data, coord, cov.mod, loc.model, scale.model, shape.model,
                             print.start.values = TRUE, method = "Nelder",
                             ...){

  param <- c("sill", "range", "smooth")
  fixed.param <- list(...)[names(list(...)) %in% param]
  
  if (print.start.values)
    cat("Computing appropriate starting values\n")
  
  n.site <- ncol(data)

  loc <- scale <- shape <- rep(NA, n.site)
  dataFrech <- data
  for (i in 1:n.site){
    marg.param <- gevmle(data[,i])
    loc[i] <- marg.param["loc"]
    scale[i] <- marg.param["scale"]
    shape[i] <- marg.param["shape"]
    dataFrech[,i] <- gev2frech(dataFrech[,i], loc[i], scale[i], shape[i])
  }

  locCoeff <- loc.model$init.fun(loc)
  scaleCoeff <- scale.model$init.fun(scale)
  shapeCoeff <- shape.model$init.fun(shape)
  
  locCoeff[is.na(locCoeff)] <- 0
  scaleCoeff[is.na(scaleCoeff)] <- 0
  shapeCoeff[is.na(shapeCoeff)] <- 0
  
  if (length(fixed.param) > 0){
    args <- c(list(data = data, coord = coord, cov.mod = cov.mod,
                   marge = "emp"), fixed.param)
    cov.param <- do.call("fitcovariance", args)$param
  }
  
  else
    cov.param <- fitcovariance(data, coord, cov.mod, marge = "emp")$param

  cov.param <- schlatherfull(dataFrech, coord, start = as.list(cov.param),
                             cov.mod = cov.mod, fit.marge = FALSE, method = method,
                             warn = FALSE)$param
  
  start <- as.list(cov.param)

  names(locCoeff) <- names(scaleCoeff) <- names(shapeCoeff) <- NULL
  
  start <- c(start, as.list(unlist(list(locCoeff = locCoeff,
                                        scaleCoeff = scaleCoeff,
                                        shapeCoeff = shapeCoeff))))

  if (print.start.values)
    cat("Starting values are defined\n")

  if (print.start.values){
    cat("Starting values are:\n")
    print(unlist(start))
  }
  
  return(start)
}


.start.schlatherind <- function(data, coord, cov.mod, loc.model, scale.model, shape.model,
                                print.start.values = TRUE, method = "Nelder",
                                ...){

  param <- c("alpha", "sill", "range", "smooth")
  fixed.param <- list(...)[names(list(...)) %in% param]
  
  if (print.start.values)
    cat("Computing appropriate starting values\n")
  
  n.site <- ncol(data)

  loc <- scale <- shape <- rep(NA, n.site)
  dataFrech <- data
  for (i in 1:n.site){
    marg.param <- gevmle(data[,i])
    loc[i] <- marg.param["loc"]
    scale[i] <- marg.param["scale"]
    shape[i] <- marg.param["shape"]
    dataFrech[,i] <- gev2frech(dataFrech[,i], loc[i], scale[i], shape[i])
  }

  locCoeff <- loc.model$init.fun(loc)
  scaleCoeff <- scale.model$init.fun(scale)
  shapeCoeff <- shape.model$init.fun(shape)
  
  locCoeff[is.na(locCoeff)] <- 0
  scaleCoeff[is.na(scaleCoeff)] <- 0
  shapeCoeff[is.na(shapeCoeff)] <- 0
  
  if (length(fixed.param) > 0){
    args <- c(list(data = data, coord = coord, cov.mod = cov.mod,
                   marge = "emp"), fixed.param)
    cov.param <- do.call("fitcovariance", args)$param
  }
  
  else
    cov.param <- fitcovariance(data, coord, cov.mod, marge = "emp")$param

  cov.param <- schlatherindfull(dataFrech, coord, start = c(list(alpha = .5), as.list(cov.param)),
                                cov.mod = cov.mod, fit.marge = FALSE, method = method,
                                warn = FALSE)$param
  
  start <- as.list(cov.param)

  names(locCoeff) <- names(scaleCoeff) <- names(shapeCoeff) <- NULL
  
  start <- c(start, as.list(unlist(list(locCoeff = locCoeff,
                                        scaleCoeff = scaleCoeff,
                                        shapeCoeff = shapeCoeff))))

  if (print.start.values)
    cat("Starting values are defined\n")

  if (print.start.values){
    cat("Starting values are:\n")
    print(unlist(start))
  }
  
  return(start)
}


.start.schlatherind <- function(data, coord, cov.mod, loc.model, scale.model, shape.model,
                                print.start.values = TRUE, method = "Nelder",
                                ...){

  param <- c("alpha", "sill", "range", "smooth")
  fixed.param <- list(...)[names(list(...)) %in% param]
  
  if (print.start.values)
    cat("Computing appropriate starting values\n")
  
  n.site <- ncol(data)

  loc <- scale <- shape <- rep(NA, n.site)
  dataFrech <- data
  for (i in 1:n.site){
    marg.param <- gevmle(data[,i])
    loc[i] <- marg.param["loc"]
    scale[i] <- marg.param["scale"]
    shape[i] <- marg.param["shape"]
    dataFrech[,i] <- gev2frech(dataFrech[,i], loc[i], scale[i], shape[i])
  }

  locCoeff <- loc.model$init.fun(loc)
  scaleCoeff <- scale.model$init.fun(scale)
  shapeCoeff <- shape.model$init.fun(shape)
  
  locCoeff[is.na(locCoeff)] <- 0
  scaleCoeff[is.na(scaleCoeff)] <- 0
  shapeCoeff[is.na(shapeCoeff)] <- 0
  
  if (length(fixed.param) > 0){
    args <- c(list(data = data, coord = coord, cov.mod = cov.mod,
                   marge = "emp"), fixed.param)
    cov.param <- do.call("fitcovariance", args)$param
  }
  
  else
    cov.param <- fitcovariance(data, coord, cov.mod, marge = "emp")$param

  cov.param <- geomgaussfull(dataFrech, coord, start = c(list(alpha = .5), as.list(cov.param)),
                             cov.mod = cov.mod, fit.marge = FALSE, method = method,
                             warn = FALSE)$param
  
  start <- as.list(cov.param)

  names(locCoeff) <- names(scaleCoeff) <- names(shapeCoeff) <- NULL
  
  start <- c(start, as.list(unlist(list(locCoeff = locCoeff,
                                        scaleCoeff = scaleCoeff,
                                        shapeCoeff = shapeCoeff))))

  if (print.start.values)
    cat("Starting values are defined\n")

  if (print.start.values){
    cat("Starting values are:\n")
    print(unlist(start))
  }
  
  return(start)
}

.start.geomgauss <- function(data, coord, cov.mod, loc.model, scale.model, shape.model,
                             print.start.values = TRUE, method = "Nelder",
                             ...){

  param <- c("sigma2", "sill", "range", "smooth")
  fixed.param <- list(...)[names(list(...)) %in% param]
  
  if (print.start.values)
    cat("Computing appropriate starting values\n")
  
  n.site <- ncol(data)

  loc <- scale <- shape <- rep(NA, n.site)
  dataFrech <- data
  for (i in 1:n.site){
    marg.param <- gevmle(data[,i])
    loc[i] <- marg.param["loc"]
    scale[i] <- marg.param["scale"]
    shape[i] <- marg.param["shape"]
    dataFrech[,i] <- gev2frech(dataFrech[,i], loc[i], scale[i], shape[i])
  }

  locCoeff <- loc.model$init.fun(loc)
  scaleCoeff <- scale.model$init.fun(scale)
  shapeCoeff <- shape.model$init.fun(shape)
  
  locCoeff[is.na(locCoeff)] <- 0
  scaleCoeff[is.na(scaleCoeff)] <- 0
  shapeCoeff[is.na(shapeCoeff)] <- 0
  
  if (length(fixed.param) > 0){
    args <- c(list(data = data, coord = coord, cov.mod = cov.mod,
                   marge = "emp"), fixed.param)
    cov.param <- do.call("fitcovariance", args)$param
  }
  
  else
    cov.param <- fitcovariance(data, coord, cov.mod, marge = "emp")$param

  cov.param <- geomgaussfull(dataFrech, coord, start = c(list(sigma2 = 1), as.list(cov.param)),
                             cov.mod = cov.mod, fit.marge = FALSE, method = method,
                             warn = FALSE)$param
  
  start <- as.list(cov.param)

  names(locCoeff) <- names(scaleCoeff) <- names(shapeCoeff) <- NULL
  
  start <- c(start, as.list(unlist(list(locCoeff = locCoeff,
                                        scaleCoeff = scaleCoeff,
                                        shapeCoeff = shapeCoeff))))

  if (print.start.values)
    cat("Starting values are defined\n")

  if (print.start.values){
    cat("Starting values are:\n")
    print(unlist(start))
  }
  
  return(start)
}
