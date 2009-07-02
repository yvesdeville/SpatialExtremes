print.spatgev <- function(x, digits = max(3, getOption("digits") - 3), ...){

  cat("   Deviance:", x$deviance, "\n")

  param <- x$fitted.values
  loc.idx <- which(substr(names(param), 1, 3) == "loc")
  scale.idx <- which(substr(names(param), 1, 5) == "scale")
  shape.idx <- which(substr(names(param), 1, 5) == "shape")

  cat("    Location Parameters:\n")
  print.default(format(param[loc.idx], digits = digits), print.gap = 2, 
                quote = FALSE)
  cat("       Scale Parameters:\n")
  print.default(format(param[scale.idx], digits = digits), print.gap = 2, 
                quote = FALSE)
  cat("       Shape Parameters:\n")
  print.default(format(param[shape.idx], digits = digits), print.gap = 2, 
                quote = FALSE)
    
  if(!is.null(x$std.err)) {
    cat("\nStandard Errors\n")
    print.default(format(x$std.err, digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  if(!is.null(x$var.cov)) {
    cat("\nAsymptotic Variance Covariance\n")
    print.default(format(x$var.cov, digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  cat("\nOptimization Information\n")
  cat("  Convergence:", x$convergence, "\n")
  cat("  Function Evaluations:", x$counts["function"], "\n")
  if(!is.na(x$counts["gradient"]))
    cat("  Gradient Evaluations:", x$counts["gradient"], "\n")
  if(!is.null(x$message)) cat("\nMessage:", x$message, "\n")
  cat("\n")
}

print.maxstab <- function(x, digits = max(3, getOption("digits") - 3), ...){

  cat("        Estimator:", x$est, "\n")
  cat("            Model:", x$model, "\n")
  if (x$est == 'MPLE'){
    cat("   Pair. Deviance:", x$deviance, "\n")
    cat("              TIC:", TIC(x), "\n")
  }
  if (x$est == "Least Square")
    cat("  Objective Value:", x$opt.value, "\n")
  if ((x$model == "Schlather") || (x$model == "Geometric")){

    if (x$cov.mod == "emp")
      cov.mod <- "Empirical"
    
    if (x$cov.mod == "whitmat")
      cov.mod <- "Whittle-Matern"
    
    if (x$cov.mod == "powexp")
      cov.mod <- "Powered Exponential"

    if (x$cov.mod == "cauchy")
      cov.mod <- "Cauchy"
    
    
    cat("Covariance Family:", cov.mod, "\n")
    
    cat("\nEstimates\n")
    cat("  Marginal Parameters:\n")

    if (x$fit.marge){
      idx <- which(names(x$fitted.values) == "alpha")
      idx <- c(idx, which(names(x$fitted.values) == "sigma2"))
      idx <- c(idx, which(names(x$fitted.values) == "sill"))
      idx <- c(idx, which(names(x$fitted.values) == "range"))
      idx <- c(idx, which(names(x$fitted.values) == "smooth"))
      

      margin.param <- x$fitted.values[-idx]
      loc.idx <- which(substr(names(margin.param), 1, 3) == "loc")
      scale.idx <- which(substr(names(margin.param), 1, 5) == "scale")
      shape.idx <- which(substr(names(margin.param), 1, 5) == "shape")

      cat("    Location Parameters:\n")
      print.default(format(margin.param[loc.idx], digits = digits), print.gap = 2, 
                    quote = FALSE)
      cat("       Scale Parameters:\n")
      print.default(format(margin.param[scale.idx], digits = digits), print.gap = 2, 
                    quote = FALSE)
      cat("       Shape Parameters:\n")
      print.default(format(margin.param[shape.idx], digits = digits), print.gap = 2, 
                    quote = FALSE)
      cat("  Dependence Parameters:\n")
      print.default(format(x$fitted.values[idx], digits = digits), print.gap = 2, 
                    quote = FALSE)
    }

    else{
      cat("  Assuming unit Frechet.\n")
      cat("  Dependence Parameters:\n")
      print.default(format(x$fitted.values, digits = digits), print.gap = 2, 
                    quote = FALSE)
    }
  }

  else{
    cat("Covariance Family:", x$cov.mod, "\n")
    
    cat("\nEstimates\n")
    cat("  Marginal Parameters:\n")

    if (x$fit.marge){
      idx <- which(substr(names(x$fitted.values), 1, 3) == "cov")

      margin.param <- x$fitted.values[-idx]
      loc.idx <- which(substr(names(margin.param), 1, 3) == "loc")
      scale.idx <- which(substr(names(margin.param), 1, 5) == "scale")
      shape.idx <- which(substr(names(margin.param), 1, 5) == "shape")

      cat("    Location Parameters:\n")
      print.default(format(margin.param[loc.idx], digits = digits), print.gap = 2, 
                    quote = FALSE)
      cat("       Scale Parameters:\n")
      print.default(format(margin.param[scale.idx], digits = digits), print.gap = 2, 
                    quote = FALSE)
      cat("       Shape Parameters:\n")
      print.default(format(margin.param[shape.idx], digits = digits), print.gap = 2, 
                    quote = FALSE)
      cat("  Dependence Parameters:\n")
      print.default(format(x$fitted.values[idx], digits = digits), print.gap = 2, 
                    quote = FALSE)
    }

    else{
      cat("  Not estimated.\n")
      cat("  Dependence Parameters:\n")
      print.default(format(x$fitted.values, digits = digits), print.gap = 2, 
                    quote = FALSE)
    }
  }
    
  if(!is.null(x$std.err)) {
    cat("\nStandard Error Type:", x$std.err.type, "\n")
    cat("\nStandard Errors\n")
    print.default(format(x$std.err, digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  if(!is.null(x$var.cov)) {
    cat("\nAsymptotic Variance Covariance\n")
    print.default(format(x$var.cov, digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  if(!is.null(x$corr)) {
    cat("\nCorrelation\n")
    print.default(format(x$corr, digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  cat("\nOptimization Information\n")
  cat("  Convergence:", x$convergence, "\n")
  cat("  Function Evaluations:", x$counts["function"], "\n")
  if(!is.na(x$counts["gradient"]))
    cat("  Gradient Evaluations:", x$counts["gradient"], "\n")
  if(!is.null(x$message)) cat("\nMessage:", x$message, "\n")
  cat("\n")
}

logLik.maxstab <- function(object, ...){
  llk <- object$logLik
  attr(llk, "df") <- length(fitted(object))
  class(llk) <- "logLik"
  return(llk)
}

profile2d <- function(fitted, ...){
  UseMethod("profile2d")
}

print.pspline <- function(x, ...){
  cat("Call:\n")
  print(x$call)

  cat("\n  Rank:", x$rank, "\t(G)CV Score:", round(x$cv, 3),
      "\n")
  cat("Degree:", x$degree, "\t Penalty: ",
      round(x$penalty, 3), "\n")
  cat("\n     Degree of freedom:", round(x$df, 3), "\n")
  cat("Res. Degree of freedom:", round(x$res.df, 3), "\n")  
}

TIC <- function(object, ...){
  UseMethod("TIC")
}
