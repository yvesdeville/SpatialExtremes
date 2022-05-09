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

## *****************************************************************************
## Overload the existing predict method. 
## XXX TODO
## arrange the output. Remove the possibility of several confidence
## levels and use a dimension c("L" , "mean" , "U") instead of two.
## So the output will be a four-dimensional array with dimensions
## site, obs, output type and return period.
##
## *****************************************************************************

##'
##' @description Predict method for the \code{"spatgev"} class.
##'
##' @details The output depends on the choice given by the formal
##'     argument \code{out}. In both cases this is a list with two
##'     elements \code{GEVParam} and \code{RL} giving the GEV
##'     parameters and the Return Level (provided that \code{ret.lev}
##'     was not of length zero.
##'
##' \itemize{ 
##'     \item{\code{"data"} }{
##'         The elements \code{GEVParam} and \code{RL} ate \emph{data
##'         frames}. Both have columns with understandable names such
##'         as \code{"Quant"} for the return level a.k.a. quantile.
##'         This choice is relevant when plots are to be produced
##'         using \pkg{ggplot2}. The variables specified with
##'         \code{keep} can then be used to define the aesthetics or
##'         the facets in trellis graphics.
##'      }
##'      \item{\code{"array"} }{
##'         The elements \code{GEVParam} and \code{RL} arre two
##'         \emph{arrays}. In both cases the first two dimensions
##'         match the spatial and the temporal dimensions. The other
##'         dimensions give the the return period, the return level,
##'         the type of return level (estimate, lower upper bound) or
##'         the GEV parameter (with value \code{"loc"} \code{"scale"}
##'         or \code{"shape"}).
##'     }
##' }
##' 
##' @title Predict Method for the \code{"spatgev"} Class
##' 
##' @param object An object with class \code{"spatgev"}.
##' 
##' @param newdata A matrix or data frame with the new "spatial"
##'     covariates.
##'
##' @param temp.newdata An optional matrix or data frame with the new
##'     "temporal" covariates.
##'
##' @param keep,temp.keep Optional character vectors giving the name
##'     of columns to be kept in the output. These will actually be
##'     used only when \code{out} specifies a data frame output
##'     i.e. when \code{out} is \code{"data"}. \bold{Not implemented
##'     yet}.
##' 
##' @param ret.per A vector of return periods.
##' 
##' @param level A confidence level between 0.0 and 1.0. If
##'     \code{NULL} of of length zero the confidence limits will not
##'     be computed.
##' 
##' @param out Character telling what kind of output is wanted. See
##'     \bold{Details}.
##'
##' @param ... Not used yet.
##'
##' @return A data frame or a list of arrays depending on the value of
##'     \code{out}, see \bold{Details}.
##'
##' @section Caution: This method will work correctly only it at the
##'     creation of \code{object} the values 
##'     \itemize{
##'     \item{\code{data} }{
##'         Should have rownames allowing the identification of the
##'          observations. These can be years or dates in POSIX format
##'          like \code{"2020-01-01"}. Should also have colnames
##'          allowing the identification of the sites.
##'     }
##'     \item{\code{covariables} }{
##'         Should have rownames allowing the identification of the
##'         sites, identical to the colnames of \code{data}.
##'     }
##'    \item{\code{temp.data} } {
##'        If used whould have rownames allowing the identification of
##'      the observations, identical to the rownames of \code{data}.
##'    }
##'    }
##'    The same rules should apply to \code{newdata} and \code{temp.newdata}.
##'    Note that when \code{temp.data} was used at the creation of \code{object}
##'    \code{temp.newdata} should be given.
##' 
predict.spatgev <- function(object,
                            newdata,
                            temp.newdata,
                            keep = NULL,
                            temp.keep = NULL,
                            ret.per = c(10, 50, 100),
                            level = 0.95,
                            out = c("data", "array"),
                            ...) {

    GEVNames <- c("loc", "scale", "shape")
    param <- object$param
    out <- match.arg(out)
    
    ## data
    if (!missing(newdata)){
        spatData <- newdata
        if (is.null(dim(spatData))) spatData <- t(as.matrix(spatData))
    } else spatData <- object$covariables
    
    ## =========================================================================
    ## Matrix of "temporal covariates data". If no temporal covariate
    ## is used we still define a matrix with one row in order to avoid
    ## multiple exceptions later. So in this case 'tempnNew' is set to
    ## 1. If no temporal covariate is used we may consider that the
    ## prediction is still for a specific time block.
    ## =========================================================================
    
    if (any(object$use.temp.cov)) {
        if (!missing(temp.newdata)) {
            tempData <- temp.newdata
            if (is.null(dim(tempData))) tempData <- t(as.matrix(tempData))
        } else {
            ## temp.data <- object$temp.cov[nrow(object$temp.cov), , drop = FALSE]
            tempData <- object$temp.cov
        }
        ## tempnNew <- nrow(tempData)
    } else {
        tempData <- matrix(as.numeric(NA), nrow = 1, ncol = 1)
    }
    spatnNew <- nrow(spatData)
    tempnNew <- nrow(tempData)
    
    ## =========================================================================
    ## Find the indices from the prefixes
    ## =========================================================================
    
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
    ## Store the design matrices in two lists for later use, and
    ## compute the GEV parameters 'theta' or the gradients:
    ##
    ## o The 'spatX' design matrix has 'spatnNew' rows,
    ##
    ## o The 'tempX' design matrix has 'tempnNew' rows.
    ## 
    ## =========================================================================

    spatX <- tempX <- list()
    theta <- array(0.0, dim = c(spatnNew, tempnNew, 3),
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
      
        ## =====================================================================
        ## Compute the conficence intervals on the return levels.
        ## =====================================================================
        
        if (length(level)) {
            
            if (length(level) > 1) stop("'level' must be of length 0 or 1")
            fLevel <- as.character(level)
            RL <- array(0.0,
                        dim = c(spatnNew, tempnNew, length(ret.per), 3L),
                        dimnames = list("site" = rownames(spatData),
                                        "obs" = rownames(tempData),
                                        "period" = as.character(ret.per),
                                        "which" = c("Est", "L", "U")))
            probL <- (1 - level) / 2
            probU <- 1 - probL
            q <- qnorm(c("L" = probL, "U" = probU), mean = 0.0, sd = 1.0)
            
            for (ell in seq_along(ret.per)) {
                
                vec <- .qgevDer(1 - 1 / ret.per[ell],
                                theta[ , , "loc"],
                                theta[ , , "scale"],
                                theta[ , , "shape"])
                RL[ , , ell, "Est"] <- vec
                
                ## =============================================================
                ## 'Jtheta' and 'Jpsi' are arrays of derivatives of the
                ## return levels at the chosen sites ans obs.  Their first
                ## two dimensions being 'site' and 'obs'.
                ##
                ## o 'Jtheta' contains the derivatives w.r.t. the
                ## GEV parameters 'theta',
                ## 
                ## o 'Jpsi' contains the derivatives w.r.t the p
                ## parameters of the model.
                ## =============================================================
                
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
                        RL[i, j, ell, c("L", "U")] <-  RL[i, j, ell, "Est"] + q * sqrt(v)
                    }
                }
                
            }
            
        } else {
            RL <- array(0.0, dim = c(spatnNew, tempnNew, length(ret.per)),
                        dimnames = list(rownames(spatData), rownames(tempData),
                                        as.character(ret.per)))
            for (ell in seq_along(ret.per)) {
                RL[ , , ell] <- .qgev(1 - 1 / ret.per[ell], theta[ , , "loc"],
                                      theta[ , , "scale"], theta[ , , "shape"])
            }   
        }
    } else RL <- NULL

        
    if (out == "array") {
        return(list(GEVParam = theta, RL = RL))
    } else {

        ### Parameters
        GEVParam <- aperm(theta, perm = c(2, 1, 3))
        dim(GEVParam) <- c(spatnNew * tempnNew, 3)
        colnames(GEVParam) <- GEVNames
        
        if (is.null(rownames(spatData))) {
            Where <- rep(as.character(NA), nrow(GEVParam))
        } else {
            Where <- rep(rownames(spatData), each = nrow(tempData))
        }
        if (is.null(rownames(tempData))) {
            When <- rep(as.character(NA), nrow(GEVParam))
        } else {
            When <- rep(rownames(tempData), time = nrow(spatData))
        }
        
        GEVParam <- data.frame(Where, When, GEVParam)

        ## Return levels
        if (length(ret.per)) {

            if (length(level)) {
                RL <- aperm(RL, perm = c(2, 1, 3, 4))
                RLnms <- paste("RL", dimnames(RL)[[4]], sep = "_")
                dim(RL) <- c(spatnNew * tempnNew, dim(RL)[3], dim(RL)[4])
            } else {
                RL <- aperm(RL, perm = c(2, 1, 3))
                dim(RL) <- c(spatnNew * tempnNew, dim(RL)[3], 1)
                RLnms <- "RL_Est"
            }
            RLOut <- data.frame(Where, When, Period = ret.per[1], RL[ , 1, ])
            
            colnames(RLOut) <- c("Where", "When", "Period", RLnms)
            
            for (ell in 2:length(ret.per)) {
                df <- data.frame(Where, When, Period = ret.per[ell], RL[ , ell, ])
                colnames(df) <- c("Where", "When", "Period", RLnms)
                RLOut <- rbind(RLOut, df)
            }
        } else {
            RLOut <- NULL
        }
             
        list(GEVParam = GEVParam, RL = RLOut) 
        
    }
        


        
}
