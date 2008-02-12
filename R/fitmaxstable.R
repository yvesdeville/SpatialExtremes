fitmaxstab <- function(data, coord, cov.mod, loc.form, scale.form, shape.form,
                       marg.cov = NULL, iso = FALSE, ..., fit.marge = FALSE, warn = TRUE,
                       method = "Nelder", start, control = list(),
                       std.err.type = "score", corr = FALSE){

  if (!(std.err.type) %in% c("none", "score", "grad"))
    stop("'std.err.type' must be one of 'none', 'score' or 'grad'")
  
  if (nrow(coord) != ncol(data))
    stop("'data' and 'coord' don't match")

  if (!is.null(marg.cov) && is.null(colnames(marg.cov)))
    stop("'marg.cov' must have named columns")

  if (!is.null(marg.cov) && (nrow(marg.cov) != nrow(coord)))
    stop("'data' and 'marg.cov' don't match")
      
  if (missing(loc.form) && missing(scale.form) && missing(shape.form))
    reg.mod <- "full"
  
  if (!missing(loc.form) && !missing(scale.form) && !missing(shape.form)){
    reg.mod <- "spatgev"
    fit.marge <- TRUE
    
    if ((class(loc.form) != "formula") || (class(scale.form) != "formula") ||
        (class(shape.form) != "formula"))
      stop("''loc.form'', ''scale.form'' and ''shape.form'' must be valid R formulas")
  }
  
  flag <- missing(loc.form) + missing(scale.form)  + missing(shape.form)
  
  if (!(flag %in% c(0, 3)))
    stop("if one formula is given for the GEV parameters, then it should
be given for *ALL* GEV parameters")

  if (method != "nlminb"){
    if (is.null(control$maxit))
      control$maxit <- 10000
  }

  else{
    if (is.null(control$eval.max))
      control$eval.max <- 15000
    
    if (is.null(control$iter.max))
      control$iter.max <- 10000
  }
  
  if (cov.mod == "gauss")
    fitted <- switch(reg.mod, "full" = smithfull(data, coord, ..., fit.marge = fit.marge,
                                iso = iso, warn = warn, method = method, control = control,
                                std.err.type = std.err.type, corr = corr, start = start),
                     "spatgev" = smithform(data, coord, ..., loc.form = loc.form, scale.form = scale.form,
                       shape.form = shape.form, fit.marge = fit.marge, iso = iso, marg.cov = marg.cov,
                       warn = warn, method = method, control = control, std.err.type =
                       std.err.type, corr = corr, start = start))
  
  
  else{

    if (substr(cov.mod, 1, 1) == "i")
      fitted <- switch(reg.mod, "full" = schlatherindfull(data, coord, cov.mod = substr(cov.mod, 2, 8),
                                    ..., fit.marge = fit.marge, warn = warn,
                                    method = method, control = control, std.err.type = std.err.type,
                                    corr = corr, start = start),
                         "spatgev" = schlatherindform(data, coord, cov.mod = substr(cov.mod, 2, 8), ...,
                           loc.form = loc.form, scale.form = scale.form, shape.form = shape.form,
                           fit.marge = fit.marge, marg.cov = marg.cov, warn = warn,
                           method = method, control = control, std.err.type = std.err.type, corr = corr,
                           start = start))

      else {
        if (substr(cov.mod, 1, 1) == "g")
          fitted <- switch(reg.mod, "full" = geomgaussfull(data, coord, cov.mod = substr(cov.mod, 2, 8),
                                      ..., fit.marge = fit.marge, warn = warn,
                                      method = method, control = control, std.err.type = std.err.type,
                                      corr = corr, start = start),
                           "spatgev" = geomgaussform(data, coord, cov.mod = substr(cov.mod, 2, 8), ...,
                             loc.form = loc.form, scale.form = scale.form, shape.form = shape.form,
                             fit.marge = fit.marge, marg.cov = marg.cov, warn = warn,
                             method = method, control = control, std.err.type = std.err.type, corr = corr,
                             start = start))

        else
          fitted <- switch(reg.mod, "full" = schlatherfull(data, coord, cov.mod = cov.mod,
                                      ..., fit.marge = fit.marge, warn = warn,
                                      method = method, control = control, std.err.type = std.err.type,
                                      corr = corr, start = start),
                           "spatgev" = schlatherform(data, coord, cov.mod = cov.mod, ...,
                             loc.form = loc.form, scale.form = scale.form, shape.form = shape.form,
                             fit.marge = fit.marge, marg.cov = marg.cov, warn = warn,
                             method = method, control = control, std.err.type = std.err.type, corr = corr,
                             start = start))
      }
  }
  
  return(fitted)
}

fitnsmaxstab <- function(data, coord, cov.mod, sigma2.form, loc.form, scale.form, shape.form,
                         marg.cov = NULL, ..., fit.marge = FALSE, warn = TRUE,
                         method = "Nelder", start, control = list(),
                         std.err.type = "score", corr = FALSE){

  if (!(std.err.type) %in% c("none", "score", "grad"))
    stop("'std.err.type' must be one of 'none', 'score' or 'grad'")
  
  if (nrow(coord) != ncol(data))
    stop("'data' and 'coord' don't match")

  if (!is.null(marg.cov) && is.null(colnames(marg.cov)))
    stop("'marg.cov' must have named columns")

  if (!is.null(marg.cov) && (nrow(marg.cov) != nrow(coord)))
    stop("'data' and 'marg.cov' don't match")
      
  if (missing(loc.form) && missing(scale.form) && missing(shape.form))
    reg.mod <- "full"
  
  if (!missing(loc.form) && !missing(scale.form) && !missing(shape.form)){
    reg.mod <- "spatgev"
    fit.marge <- TRUE
    
    if ((class(loc.form) != "formula") || (class(scale.form) != "formula") ||
        (class(shape.form) != "formula"))
      stop("''loc.form'', ''scale.form'' and ''shape.form'' must be valid R formulas")
  }
  
  flag <- missing(loc.form) + missing(scale.form)  + missing(shape.form)
  
  if (!(flag %in% c(0, 3)))
    stop("if one formula is given for the GEV parameters, then it should
be given for *ALL* GEV parameters")

  if (method != "nlminb"){
    if (is.null(control$maxit))
      control$maxit <- 10000
  }

  else{
    if (is.null(control$eval.max))
      control$eval.max <- 15000
    
    if (is.null(control$iter.max))
      control$iter.max <- 10000
  }
  
  fitted <- switch(reg.mod, "full" = nsgeomgaussfull(data, coord, cov.mod = cov.mod,
                              sigma2.form, ..., fit.marge = fit.marge, warn = warn,
                              method = method, control = control, std.err.type = std.err.type,
                              corr = corr, start = start),
                   "spatgev" = nsgeomgaussform(data, coord, cov.mod = cov.mod,
                     sigma2.form, ..., loc.form = loc.form, scale.form = scale.form,
                     shape.form = shape.form, fit.marge = fit.marge, marg.cov = marg.cov, warn = warn,
                     method = method, control = control, std.err.type = std.err.type, corr = corr,
                     start = start))

  return(fitted)
}
