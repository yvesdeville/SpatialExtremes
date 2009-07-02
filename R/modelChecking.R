qqextcoeff <- function(fitted, estim = "ST", marge = "emp",
                       xlab = "Semi-Empirical", ylab = "Model", ...){

  if (!any("maxstab" %in% class(fitted)))
    stop("This functin is only available for 'maxstab' objects")

  model <- fitted$model
  data <- fitted$data
  coord <- fitted$coord
  ext.coeff <- fitted$ext.coeff

  if (model == "Smith"){
    dist <- distance(coord, vec = TRUE)
    exco.mod <- apply(dist, 1, ext.coeff)
  }

  else{
    dist <- distance(coord)
    exco.mod <- ext.coeff(dist)
  }

  exco.emp <- fitextcoeff(data, coord, plot = FALSE, estim = estim,
                          marge = marge)$ext.coeff[,"ext.coeff"]

  plot(exco.emp, exco.mod, xlab = xlab, ylab = ylab, ...)
  abline(0, 1)
}

qqgev <- function(fitted, xlab, ylab, ...){

  data <- fitted$data
  n.site <- ncol(data)

  gev.param <- matrix(NA, ncol = 3, nrow = n.site)
  colnames(gev.param) <- c("loc", "scale", "shape")

  for (i in 1:n.site)
    gev.param[i,] <- gevmle(data[,i])

  pred <- predict(fitted)

  if (missing(xlab))
    xlab <- c(expression(mu[MLE]), expression(sigma[MLE]), expression(xi[MLE]))

  if (missing(ylab))
    ylab <- c(expression(mu[Model]), expression(sigma[Model]), expression(xi[Model]))

  op <- par(mfrow=c(1,3))
  on.exit(par(op))

  if (length(unique(pred[,"loc"])) != 1){
    plot(gev.param[,"loc"], pred[,"loc"], xlab = xlab[1], ylab = ylab[1], ...)
    abline(0, 1)
  }

  else{
    hist(gev.param[,"loc"], xlab = xlab[1], main = "")
    axis(3, at = pred[1, "loc"], labels = ylab[1])
  }

  if (length(unique(pred[,"scale"])) !=1){
    plot(gev.param[,"scale"], pred[,"scale"], xlab = xlab[2], ylab = ylab[2], ...)
    abline(0, 1)
  }

  else{
    hist(gev.param[,"scale"], xlab = xlab[2], main = "")
    axis(3, at = pred[1, "scale"], labels = ylab[2])
  }

  if (length(unique(pred[,"shape"])) != 1){
    plot(gev.param[,"shape"], pred[,"shape"], xlab = xlab[3], ylab = ylab[3], ...)
    abline(0, 1)
  }

  else{
    hist(gev.param[,"shape"], xlab = xlab[3], main = "")
    axis(3, at = pred[1, "shape"], labels = ylab[3])
  }
}
