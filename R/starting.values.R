.start.smith <- function(data, coord, covariables, loc.form, scale.form, shape.form,
                         print.start.values = TRUE, method = "Nelder",
                         iso = TRUE, ...){

  n.site <- ncol(data)
  
  if (ncol(coord) == 2)
    param <- c("cov11", "cov12", "cov22")

  else
    param <- c("cov11", "cov12", "cov13", "cov22", "cov23", "cov33")

  fixed.param <- list(...)[names(list(...)) %in% param]
  idx.cov <- which(names(fixed.param) %in% param)
  fixed.param.cov <- fixed.param[idx.cov]
  fixed.param.gev <- fixed.param[-idx.cov]
  
  if (print.start.values)
    cat("Computing appropriate starting values\n")
  
  if (length(fixed.param.cov) > 0){
    args <- c(list(data = data, coord = coord, marge = "emp", iso = iso),
              fixed.param.cov)
    covs <- do.call("fitcovmat", args)$fitted
  }

  else
    covs <- fitcovmat(data, coord, marge = "emp", iso = iso)$fitted

  if (iso){
    covs <- covs[1]
    names(covs) <- "cov"
  }

  args <- c(list(data = data, covariables = as.matrix(covariables), loc.form = loc.form,
                 scale.form = scale.form, shape.form = shape.form, std.err.type = "none",
                 method = method), fixed.param.gev)

  spatgev <- do.call("fitspatgev", args)

  frech <- data
  gev <- predict(spatgev)
  for (i in 1:n.site)
    frech[,i] <- gev2frech(frech[,i], gev[i,"loc"], gev[i,"scale"], gev[i,"shape"])

  args <- c(list(data = frech, coord = coord, fit.marge = FALSE, start = as.list(covs),
                 iso = iso, warn = FALSE, method = method, std.err.type = "none"),
            fixed.param.cov)

  covs <- do.call("smithfull", args)
    
  start <- c(as.list(covs$param), as.list(spatgev$param))

  if (print.start.values){
    cat("Starting values are defined\n")
    cat("Starting values are:\n")
    print(c(covs$fitted, spatgev$fitted))
  }

  return(start)
}


.start.schlather <- function(data, coord, covariables, cov.mod, loc.form,
                             scale.form, shape.form, print.start.values = TRUE,
                             method = "Nelder", ...){

  n.site <- ncol(data)
  param <- c("sill", "range", "smooth")
  fixed.param <- list(...)[names(list(...)) %in% param]

  idx.cov <- which(names(fixed.param) %in% param)
  fixed.param.cov <- fixed.param[idx.cov]
  fixed.param.gev <- fixed.param[-idx.cov]
  
  if (print.start.values)
    cat("Computing appropriate starting values\n")
  
  if (length(fixed.param.cov) > 0){
    args <- c(list(data = data, coord = coord, cov.mod = cov.mod,
                   marge = "emp"), fixed.param.cov)
    cov.param <- do.call("fitcovariance", args)$fitted
  }
  
  else
    cov.param <- fitcovariance(data, coord, cov.mod, marge = "emp")$fitted

  args <- c(list(data = data, covariables = as.matrix(covariables), loc.form = loc.form,
                 scale.form = scale.form, shape.form = shape.form, std.err.type = "none",
                 method = method), fixed.param.gev)

  spatgev <- do.call("fitspatgev", args)

  frech <- data
  gev <- predict(spatgev)
  for (i in 1:n.site)
    frech[,i] <- gev2frech(frech[,i], gev[i,"loc"], gev[i,"scale"], gev[i,"shape"])

  args <- c(list(data = frech, coord = coord, cov.mod = cov.mod, fit.marge = FALSE,
                 start = as.list(cov.param), warn = FALSE, method = method,
                 std.err.type = "none"), fixed.param.cov)

  cov.param <- do.call("schlatherfull", args)
  
  start <- c(as.list(cov.param$param), as.list(spatgev$param))

  if (print.start.values){
    cat("Starting values are defined\n")
    cat("Starting values are:\n")
    print(c(cov.param$fitted, spatgev$fitted))
  }
  
  return(start)
}


.start.schlatherind <- function(data, coord, covariables, cov.mod, loc.form,
                                scale.form, shape.form, print.start.values = TRUE,
                                method = "Nelder", ...){

  n.site <- ncol(data)
  param <- c("alpha", "sill", "range", "smooth")
  fixed.param <- list(...)[names(list(...)) %in% param]

  idx.cov <- which(names(fixed.param) %in% param)
  fixed.param.cov <- fixed.param[idx.cov]
  fixed.param.gev <- fixed.param[-idx.cov]
  
  if (print.start.values)
    cat("Computing appropriate starting values\n")
  
  if (length(fixed.param.cov) > 0){
    args <- c(list(data = data, coord = coord, cov.mod = paste("i", cov.mod, sep=""),
                   marge = "emp"), fixed.param.cov)
    cov.param <- do.call("fitcovariance", args)$fitted
  }
  
  else
    cov.param <- fitcovariance(data, coord, paste("i", cov.mod, sep=""),
                               marge = "emp")$fitted

  args <- c(list(data = data, covariables = as.matrix(covariables), loc.form = loc.form,
                 scale.form = scale.form, shape.form = shape.form, std.err.type = "none",
                 method = method), fixed.param.gev)

  spatgev <- do.call("fitspatgev", args)

  frech <- data
  gev <- predict(spatgev)
  for (i in 1:n.site)
    frech[,i] <- gev2frech(frech[,i], gev[i,"loc"], gev[i,"scale"], gev[i,"shape"])

  args <- c(list(data = frech, coord = coord, cov.mod = cov.mod, fit.marge = FALSE,
                 start = as.list(cov.param), warn = FALSE, method = method,
                 std.err.type = "none"), fixed.param.cov)

  cov.param <- do.call("schlatherindfull", args)
    
  start <- c(as.list(cov.param$param), as.list(spatgev$param))
  
  if (print.start.values){
    cat("Starting values are defined\n")
    cat("Starting values are:\n")
    print(c(cov.param$fitted, spatgev$fitted))
  }
  
  return(start)
}

.start.geomgauss <- function(data, coord, covariables, cov.mod, loc.form,
                             scale.form, shape.form, print.start.values = TRUE,
                             method = "Nelder", ...){

  n.site <- ncol(data)
  param <- c("sigma2", "sill", "range", "smooth")
  fixed.param <- list(...)[names(list(...)) %in% param]

  idx.cov <- which(names(fixed.param) %in% param)
  fixed.param.cov <- fixed.param[idx.cov]
  fixed.param.gev <- fixed.param[-idx.cov]
  
  if (print.start.values)
    cat("Computing appropriate starting values\n")
  
  if (length(fixed.param.cov) > 0){
    args <- c(list(data = data, coord = coord, cov.mod = paste("g", cov.mod, sep=""),
                   marge = "emp"), fixed.param.cov)
    cov.param <- do.call("fitcovariance", args)$fitted
  }
  
  else
    cov.param <- fitcovariance(data, coord, paste("g", cov.mod, sep=""),
                               marge = "emp")$fitted

  args <- c(list(data = data, covariables = as.matrix(covariables), loc.form = loc.form,
                 scale.form = scale.form, shape.form = shape.form, std.err.type = "none",
                 method = method), fixed.param.gev)

  spatgev <- do.call("fitspatgev", args)

  frech <- data
  gev <- predict(spatgev)
  for (i in 1:n.site)
    frech[,i] <- gev2frech(frech[,i], gev[i,"loc"], gev[i,"scale"], gev[i,"shape"])

  args <- c(list(data = frech, coord = coord, cov.mod = cov.mod, fit.marge = FALSE,
                 start = as.list(cov.param), warn = FALSE, method = method,
                 std.err.type = "none"), fixed.param.cov)

  cov.param <- do.call("geomgaussfull", args)
  
  start <- c(as.list(cov.param$param), as.list(spatgev$param))
  
  if (print.start.values){
    cat("Starting values are defined\n")
    cat("Starting values are:\n")
    print(c(cov.param$fitted, spatgev$fitted))
  }
  
  return(start)
}

.start.nsgeomgauss <- function(data, coord, covariables, cov.mod, loc.form,
                               scale.form, shape.form, sigma2.form,
                               print.start.values = TRUE, method = "Nelder",
                               ...){

  n.site <- ncol(data)
  sigma2.terms <- terms(sigma2.form)

  if (length(attributes(sigma2.terms)$factor) == 0)
    n.sigma2coeff <- attributes(sigma2.terms)$intercept

  else
    n.sigma2coeff <- attributes(sigma2.terms)$intercept +
      ncol(attributes(sigma2.terms)$factors)

  if (n.sigma2coeff == 1)
    sigma2.names <- "sigma2"

  else
    sigma2.names <- paste("sigma2Coeff", 1:n.sigma2coeff, sep="")

  param <- c(sigma2.names, "sill", "range", "smooth")
  fixed.param <- list(...)[names(list(...)) %in% param]

  idx.cov <- which(names(fixed.param) %in% param)
  fixed.param.cov <- fixed.param[idx.cov]
  fixed.param.gev <- fixed.param[-idx.cov]
  
  if (print.start.values)
    cat("Computing appropriate starting values\n")
  
  if (length(fixed.param.cov) > 0){
    args <- c(list(data = data, coord = coord, cov.mod = cov.mod,
                   marge = "emp"), fixed.param.cov)
    cov.param <- do.call("fitcovariance", args)$fitted
  }
  
  else
    cov.param <- fitcovariance(data, coord, cov.mod, marge = "emp")$fitted

  args <- c(list(data = data, covariables = as.matrix(covariables), loc.form = loc.form,
                 scale.form = scale.form, shape.form = shape.form, std.err.type = "none",
                 method = method), fixed.param.gev)

  spatgev <- do.call("fitspatgev", args)

  frech <- data
  gev <- predict(spatgev)
  for (i in 1:n.site)
    frech[,i] <- gev2frech(frech[,i], gev[i,"loc"], gev[i,"scale"], gev[i,"shape"])

  sigma2param <- rep(0, length(sigma2.names))
  names(sigma2param) <- sigma2.names

  args <- c(list(data = frech, coord = coord, cov.mod = cov.mod, fit.marge = FALSE,
                 sigma2.form = sigma2.form, warn = FALSE, method = method, start =
                 c(as.list(sigma2param), as.list(cov.param)), std.err.type = "none"),
            fixed.param.cov)

  cov.param <- do.call("nsgeomgaussfull", args)
  
  start <- c(as.list(cov.param$param), as.list(spatgev$param))

  if (print.start.values){
    cat("Starting values are defined\n")
    cat("Starting values are:\n")
    print(c(cov.param$fitted, spatgev$fitted))
  }
  
  return(start)
}

