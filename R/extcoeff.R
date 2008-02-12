madogram <- function(data, coord, bins, gev.param = c(0, 1, 0),
                     which = c("mado", "ext"), xlab, ylab,
                     angles = NULL, ...){
  
  if (nrow(coord) != ncol(data))
    stop("'data' and 'coord' don't match")
  
  n.site <- ncol(data)
  dist <- distance(coord)

  if (!is.null(angles)){
    distVec <- distance(coord, vec = TRUE)
    n.angles <- length(angles)
    angles.coord <- atan2(distVec[,2], distVec[,1])

    col <- rep(NA, n.site * (n.site - 1) / 2)
    idx.angles <- list()
    for (i in 2:n.angles){
      idx <- which((angles.coord < angles[i]) & (angles.coord >= angles[i-1]))
      idx.angles <- c(idx.angles, list(idx))
      col[idx] <- i-1
    }
  }

  else
    col <- 1
 
  for (i in 1:n.site){
    param <- gevmle(data[,i])
    data[,i] <- gev2frech(data[,i], param[1], param[2],
                          param[3])
    data[,i] <- frech2gev(data[,i], gev.param[1],
                          gev.param[2], gev.param[3])
  }

  if (missing(bins)){
    k <- 1
    mado <- rep(NA, length = n.site * (n.site-1)/2)
    for (i in 1:(n.site-1)){
      for (j in (i+1):n.site){
        mado[k] <- mean(abs(data[,i] - data[,j])) / 2
        k <- k + 1
      }
    }
  }

  else{
    bins <- unique(c(0, bins, Inf))  
    mado <- rep(0, length = length(bins)-1)
      
    for (k in 1:(length(bins)-1)){
      idx <- which((dist <= bins[k+1]) & (dist > bins[k]))
      
      if (length(idx)>0){
        site1 <- (idx-1) %/% n.site + 1
        site2 <- (idx-1) %% n.site + 1
        
        for (i in 1:length(idx))
          mado[k] <- mado[k] + sum(abs(data[,site1[i]] -
                                       data[,site2[i]]))
        
        mado[k] <- mado[k] / 2 / length(idx) / nrow(data)
      }

      else
        mado[k] <- NA
    }
  }

  if (gev.param[3] == 0)
    ext.coeff <- exp(mado/gev.param[2])

  else
    ext.coeff <- gev2frech(gev.param[1] + mado / gamma(1 - gev.param[3]),
                           gev.param[1], gev.param[2], gev.param[3])

  if (length(which) == 2){
    op <- par(mfrow=c(1,2))
    on.exit(par(op))
  }

  if (missing(xlab))
    xlab <- "h"

  if (missing(ylab))
      ylab <- c(expression(eta(h)), expression(theta(h)))
    
  if (any(which == "mado")){
    if (missing(bins))
      plot(dist, mado, xlab = xlab, ylab = ylab[1], col = col,
           pch = col, ...)

    else
      plot(bins[-1], mado, xlab = xlab, ylab = ylab[1], col = col,
           pch = col, ...)
  }
  
  if (any(which == "ext")){
    if (missing(bins))
      plot(dist, ext.coeff, xlab = xlab, ylab = ylab[2], col = col,
           pch = col, ...)

    else
      plot(bins[-1], ext.coeff, xlab = xlab, ylab = ylab[2], col = col,
           pch = col, ...)
  }

  if (missing(bins))
    bins <- c(0, dist)
  
  invisible(cbind(bins = bins[-1], madogram = mado, ext.coeff = ext.coeff))
}

fitextcoeff <- function(data, coord, ..., estim = "ST", marge = "emp",
                        prob = 0, plot = TRUE, loess = TRUE, method = "BFGS",
                        std.err = TRUE, xlab, ylab, angles = NULL,
                        identify = FALSE){
  
  if (nrow(coord) != ncol(data))
    stop("'data' and 'coord' don't match")
  
  if (!(estim %in% c("ST", "Smith")))
    stop("'estim' must be one of 'ST' or 'Smith'")
  
  if (!(marge %in% c("emp", "mle", "frech")))
    stop("'marge' must be one of 'emp', 'mle' or 'frech'")
  
  if ((prob < 0) || (prob >= 1))
    stop("'prob' must be in [0,1)")
  
  n.site <- ncol(data)
  n.obs <- nrow(data)
  dist <- distance(coord)
  n.pairs <- n.site * (n.site - 1) / 2

  if (!is.null(angles)){
    distVec <- distance(coord, vec = TRUE)
    n.angles <- length(angles)
    angles.coord <- atan2(distVec[,2], distVec[,1])

    col <- rep(NA, n.pairs)
    idx.angles <- list()
    for (i in 2:n.angles){
      idx <- which((angles.coord < angles[i]) & (angles.coord >= angles[i-1]))
      idx.angles <- c(idx.angles, list(idx))
      col[idx] <- i-1
    }
  }

  else
    col <- 1
  
  if (marge == "mle"){
    frech <- data
    for (i in 1:n.site){
      marg.param <- gevmle(data[,i], method = method)
      frech[,i] <- gev2frech(frech[,i], marg.param["loc"],
                             marg.param["scale"],
                             marg.param["shape"])
    }
  }
  
  if (marge == "emp") {
    frech <- data
    probs <- ppoints(n.obs, a = 0)
    for (i in 1:n.site){
      idx <- order(frech[,i])
      frech[idx,i] <- probs
    }
    
    frech <- - 1 / log(frech)
    
  }

  if (marge == "frech")
    frech <- data

  ext.coeff <- rep(NA, n.pairs)
  
  if (estim == "ST"){
    std.err <- FALSE
    if (prob == 0)
      z <- 0

    else
      z <- - 1 / log(prob)
    
    x.bar <- colMeans(1/frech)
    
    lik.fun <- function(theta)
      .C("extCoeffST", as.double(frech[,pair]), as.double(x.bar[pair]), as.double(z),
         as.double(theta), as.integer(n.obs), dns = double(1),
         PACKAGE = "SpatialExtremes")$dns
    
    k <- 1
    
    for (i in 1:(n.site - 1)){
      for (j in (i+1):n.site){
        pair <- c(i,j)
        ext.coeff[k] <- optimize(lik.fun, interval = c(1, 2))$minimum
        k <- k + 1
      }
    }
  }
  
  if (estim == "Smith") {
    ##The Smith approach needs Exponential margins instead of Frechet ones
    frech <- 1 / frech    
    
    ext.coeff <- .C("extCoeffSmith", as.double(frech), as.integer(n.obs),
                    as.integer(n.site), extCoeff = double(n.pairs),
                    PACKAGE = "SpatialExtremes")$extCoeff

    if (std.err){

      ext.coeff.std.err <- matrix(NA, ncol = n.site * (n.site - 1) / 2,
                                  nrow = n.obs)
      
      if (marge == "frech"){
        for (year in 1:n.obs){
          frech.jack <- 1 / data[-year,]
    
          ext.coeff.std.err[year,] <- .C("extCoeffSmith", as.double(frech.jack),
                                         as.integer(n.obs-1), as.integer(n.site),
                                         extCoeff = double(n.pairs),
                                         PACKAGE = "SpatialExtremes")$extCoeff
        }
      }

      if (marge == "mle"){
        for (year in 1:n.obs){
          frech.jack <- data[-year,]
          for (i in 1:n.site){
            marg.param <- gevmle(data[-year,i], method = method)
            frech.jack[,i] <- gev2frech(data[-year,i], marg.param["loc"],
                                        marg.param["scale"],
                                        marg.param["shape"])
          }

          frech.jack <- 1 / frech.jack
          
          ext.coeff.std.err[year,] <- .C("extCoeffSmith", as.double(frech.jack),
                                         as.integer(n.obs-1), as.integer(n.site),
                                         extCoeff = double(n.pairs),
                                         PACKAGE = "SpatialExtremes")$extCoeff
        }
      }

      if (marge == "emp"){
        probs <- ppoints(n.obs - 1, a = 0)
        for (year in 1:n.obs){
          frech.jack <- data[-year,]
          for (i in 1:n.site){
            idx <- order(data[-year,i])
            frech.jack[idx,i] <- probs
          }
          
          frech.jack <- - log(frech.jack)
          
          ext.coeff.std.err[year,] <- .C("extCoeffSmith", as.double(frech.jack),
                                         as.integer(n.obs-1), as.integer(n.site),
                                         extCoeff = double(n.pairs),
                                         PACKAGE = "SpatialExtremes")$extCoeff
        }
      }
      
      ext.coeff.std.err <- sqrt(apply(ext.coeff.std.err, 2, var) * (n.obs - 2) *
                                (n.obs - 1) / n.obs)
    }
  }

  if (loess){
    if (!is.null(angles)){
      loess.fit <- list()
      for (i in 1:(n.angles-1)){
        ext.coeff.sub <- ext.coeff[idx.angles[[i]]]
        dist.sub <- dist[idx.angles[[i]]]
        loess.fit <- c(loess.fit, list(loess(ext.coeff.sub ~ dist.sub)))
      }
    }

    else
      loess.fit <- loess(ext.coeff ~ dist)
  }
  
  if (plot){

    if (missing(xlab))
      xlab <- "h"

    if (missing(ylab))
      ylab <- expression(theta(h))

    if (identify)
      par(mfrow=c(1,2))
    
    plot(dist, ext.coeff, xlab = xlab, ylab = ylab, col = col, pch = col, ...)
    
    if (loess){
      h <- seq(0, max(dist), length = 200)

      if (!is.null(angles)){
        for (i in 1:(n.angles-1))
          lines(h, predict(loess.fit[[i]], data.frame(dist.sub = h)), col = i, ...)
      }

      else
        lines(h, predict(loess.fit, data.frame(dist = h)), ...)
    }

    if (identify){
      labels <- NULL
      for (i in 1:(n.site-1))
        labels <- c(labels, paste(i, "-", (i+1):n.site, sep = ""))
      
      id.points <- identify(dist, ext.coeff, labels = labels,
                            plot = FALSE)
      id.stations <- unlist(strsplit(labels[id.points], "-"))
      plot(coord)
      points(coord[as.numeric(id.stations),], col = 2)
    }
  }
  
  ans <- cbind(distance = dist, ext.coeff = ext.coeff)
  
  if (std.err)
    ans <- cbind(ans, std.err = ext.coeff.std.err)

  if (loess)
    ans <- list(ext.coeff = ans, loess = loess.fit)

  if (identify)
    ans <- c(ans, select.pairs = labels[id.points])
  
  return(ans)
}
