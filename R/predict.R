predict.maxstab <- function(object, newdata, ret.per = NULL,
                            ...){
  
  param <- object$param
    
  if (!missing(newdata)){
    data <- newdata
    marg.cov <- NULL
    
    if (is.null(dim(newdata)))
      data <- t(as.matrix(data))
  }

  else{
    data <- object$coord
    marg.cov <- object$marg.cov
  }

  if (!is.null(marg.cov))
    data <- cbind(data, marg.cov)

  if (object$fit.marge){
    loc.form <- object$loc.form
    scale.form <- object$scale.form
    shape.form <- object$shape.form
    
    loc.dsgnmat <- modeldef(data, loc.form)$dsgn.mat
    scale.dsgnmat <- modeldef(data, scale.form)$dsgn.mat
    shape.dsgnmat <- modeldef(data, shape.form)$dsgn.mat
    
    idx.loc <- which(substr(names(param), 1, 8) == "locCoeff")
    idx.scale <- which(substr(names(param), 1, 10) == "scaleCoeff")
    idx.shape <- which(substr(names(param), 1, 10) == "shapeCoeff")
    
    loc.pred <- loc.dsgnmat %*% param[idx.loc]
    scale.pred <- scale.dsgnmat %*% param[idx.scale]
    shape.pred <- shape.dsgnmat %*% param[idx.shape]
  }

  else
    loc.pred <- scale.pred <- shape.pred <- rep(1, nrow(data))

  ans <- cbind(loc.pred, scale.pred, shape.pred)
  colnames(ans) <- c("loc", "scale", "shape")
  ans <- cbind(data, ans)

  if (!is.null(ret.per)){
    ret.lev <- NULL
    for (T in ret.per)
      ret.lev <- cbind(ret.lev, .qgev(1 - 1/T, ans[,"loc"],
                                      ans[,"scale"],
                                      ans[,"shape"]))

  colnames(ret.lev) <- paste("Q", ret.per, sep="")
  ans <- cbind(ans, ret.lev)
  }
  
  return(ans)
}

predict.pspline <- function(object, new.data, ...){

  if (missing(new.data))
    new.data <- object$x

  degree <- object$degree
  knots <- object$knots
  beta <- object$beta
  
  dsgn.mat <- rb(new.data, degree = degree, knots = knots,
                 penalty = NULL)$dsgn.mat

  y <- dsgn.mat %*% beta
  ans <- cbind(new.data, y)
  colnames(ans) <- c("x", "y.hat")
  
  return(ans)
}
