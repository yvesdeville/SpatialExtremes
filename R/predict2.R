.qgevDer <- function(p, loc = 1, scale = 1, shape = 1, lower.tail = TRUE) {

    ## '.gev' will perform the error detection (p < 0, p > 1, ...)
    res <- .qgev(p, loc = loc, scale = scale, shape = shape, lower.tail = lower.tail)

    nRes <- length(res)
    ## use the same recyling as .qgev.
    p <- rep(p, length.out = nRes)
    loc <- rep(loc, length.out = nRes)
    scale <- rep(scale, length.out = nRes)
    shape <- rep(shape, length.out = nRes)
    
    Jac <- array(0.0, dim = c(length(res), 3),
                 dimnames = list(names(res), c("loc", "scale", "shape")))
    
    A <- - log(p)
    logA <- log(A)
    
    Jac[ , "loc"] <- rep(1.0, length(res))
    Jac[ , "scale"] <- -logA
    Jac[ , "shape"] <- scale * logA^2 / 2.0 

    ind <- (shape != 0.0)

    if (any(ind)) {
        Amxi <- A[ind]^(-shape[ind])
        V <- (1 - Amxi) / shape[ind]
        Jac[ind, "scale"] <- -V
        Jac[ind, "shape"] <- scale[ind] * (V - logA[ind] * Amxi) / shape[ind]
    }
       
    attr(res, "Jacobian") <- Jac
    res
    
}

## 
## XXX TODO
## arrange the output. Remove the possibility of several confidence
## levels and use a dimension c("L" , "mean" , "U") instead of two.
## So the output will be a four-dimensional array with dimensions
## site, obs, output type and return period.
## 
predict.spatgev <- function(object, newdata,
                            temp.newdata,
                            ret.per = NULL,
                            level = NULL,
                            out = c("array", "wide", "long"),
                            ...) {

    GEVNames <- c("loc", "scale", "shape")
    param <- object$param
    out <- match.arg(out)
    
    ## data
    if (!missing(newdata)){
        spatData <- newdata
        if (is.null(dim(spatData))) spatData <- t(as.matrix(spatData))
    } else spatData <- object$covariables

    if (any(object$use.temp.cov)) {
        if (!missing(temp.newdata)) {
            tempData <- temp.newdata
            if (is.null(dim(tempData))) tempData <- t(as.matrix(tempData))
        } else {
            ## temp.data <- object$temp.cov[nrow(object$temp.cov), , drop = FALSE]
            tempData <- object$temp.cov
        }
    } 

    spatnNew <- nrow(spatData)
    tempnNew <- nrow(tempData)
    
    ## Find the indices from the prefixes
    spatInd <- tempInd <- list()
    spatPrefixes <- c("loc" = "locCoeff", "scale" = "scaleCoeff",
                      "shape" = "shapeCoeff")
    tempPrefixes <- c("loc" = "tempCoeffLoc", "scale" = "tempCoeffScale",
                      "shape" = "tempCoeffShape")
    
    for (pnm in GEVNames) {
        pref <- spatPrefixes[pnm]
        idx <- which(substr(names(param), 1, nchar(pref)) == pref)
        spatInd[[pnm]] <- idx
    }
    
    for (pnm in GEVNames) {
        pref <- tempPrefixes[pnm]
        idx <- which(substr(names(param), 1, nchar(pref)) == pref)
        tempInd[[pnm]] <- idx
    }
    
    ## =========================================================================
    ## Store the designs for later use and compute the GEV parameters
    ## 'theta'.  A 'spatX' design has 'nNewSites' rows. A 'tempX'
    ## design has 'nNewObs' rows.
    ## =========================================================================
    spatX <- tempX <- list()
    theta <- array(0.0,
                   dim = c(spatnNew, tempnNew, 3),
                   dimnames = list(rownames(spatData), rownames(tempData),
                                   GEVNames))
    
    for (pnm in GEVNames) {
        spatX[[pnm]] <-
            modeldef(spatData, object[[paste0(pnm, ".form")]])$dsgn.mat
        spatTheta <- spatX[[pnm]] %*% param[spatInd[[pnm]]] 
        theta[ , , pnm] <- sweep(theta[ , , pnm, drop = FALSE],
                                 spatTheta,
                                 MARGIN = 1, FUN = "+")
        
        if (object$use.temp.cov[[pnm]]) {
            tempX[[pnm]] <- modeldef(tempData, object$temp.form[[pnm]])$dsgn.mat
            tempTheta <- tempX[[pnm]] %*% param[tempInd[[pnm]]]
            theta[ , , pnm] <-  sweep(theta[ , , pnm, drop = FALSE],
                                      tempTheta,
                                      MARGIN = 2, FUN = "+")
        } else {
            tempX[[pnm]] <- matrix(nrow = tempnNew, ncol = 0)
        }
    }
    
    ## =========================================================================
    ## Compute the return levels.
    ## =========================================================================
    
    if (!is.null(ret.per)) {
        RL <- array(0.0,
                    dim = c(spatnNew, tempnNew, length(ret.per)),
                    dimnames = list(rownames(spatData), rownames(tempData),
                                    as.character(ret.per)))
        for (ell in seq_along(ret.per)) {
            RL[ , , ell] <- .qgev(1 - 1 / ret.per[ell],
                                  theta[ , , "loc"],
                                  theta[ , , "scale"],
                                  theta[ , , "shape"])
        }
        
    }

    if (length(level)) {
        level <- sort(level)
        fLevel <- as.character(level)

        RLConf <- array(0.0,
                        dim = c(spatnNew, tempnNew, length(ret.per), 2L, length(level)),
                        dimnames = list("site" = rownames(spatData),
                                        "obs" = rownames(tempData),
                                        "period" = as.character(ret.per),
                                        "lim" = c("L", "U"),
                                        "level" = fLevel))

        probL <- (1 - level) / 2
        probU <- 1 - probL
        q <- qnorm(c("L" = probL, "U" = probU), mean = 0.0, sd = 1.0)
        
        for (k in seq_along(ret.per)) {
            
            vec <- .qgevDer(1 - 1 / ret.per[k],
                            theta[ , , "loc"],
                            theta[ , , "scale"],
                            theta[ , , "shape"])
            RL[ , , k] <- vec
        
            ## matrix with spatnNew * tempnNew rows and 3 columns
            Jtheta <- attr(vec, "Jacobian")
            dim(Jtheta) <- c(spatnNew, tempnNew, 3)
            dimnames(Jtheta) <- list(rownames(spatData), rownames(tempData),
                                  GEVNames)
            Jpsi <- array(0.0, dim = c(spatnNew, tempnNew, length(param)),
                          dimnames = list("site" = rownames(spatData),
                                          "obs" = rownames(tempData),
                                          "param" = names(param)))
                          
            for (pnm in GEVNames) {
                Jpsi[ , , spatInd[[pnm]]] <-
                    Jtheta[ , , pnm] * spatX[[pnm]][ , rep(1, tempnNew)]
                if (object$use.temp.cov[[pnm]]) {
                    Jpsi[ , , tempInd[[pnm]]] <-
                        Jtheta[ , , pnm] * t(tempX[[pnm]])[rep(1, spatnNew), ]  
                }
            }

            ## we could use some "apply" here to speed up. Worth it?
            for (i in 1:spatnNew) {
                for (j in 1:tempnNew) {
                    v <- drop(t(Jpsi[i, j, ]) %*% object$var.cov %*% Jpsi[i, j, ])
                    RLConf[i, j, k, , ell] <-  RL[i, j, k] + q * sqrt(v)
                }
            }
            
        }
    } else RLConf <- NULL
 
    return(list(RL = RL, RLConf = RLConf))
        
}
