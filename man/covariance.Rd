\name{covariance}
\alias{covariance}
\title{Defines and computes covariance functions}
\description{This function defines and computes several covariance
  function either from a fitted ``max-stable'' model; either by
  specifying directly the covariance parameters.
}
\usage{
covariance(fitted, sill, range, smooth, cov.mod = "whitmat", plot = TRUE,
dist, xlab, ylab, ...)
}
\arguments{
  \item{fitted}{An object of class ``maxstab''. Most often this will be
    the output of the \code{\link{fitmaxstab}} function. May be missing
    if \code{scale} and \code{smooth} are given.}
  \item{sill,range,smooth}{The sill, scale and smooth parameters for
    the covariance function. May be missing if \code{fitted} is given.}
  \item{cov.mod}{Character string. The name of the covariance
    model. Must be one of "whitmat", "cauchy" or "powexp" for the
    Whittle-Matern, Cauchy and Powered Exponential models.}
  \item{plot}{Logical. If \code{TRUE} (default) the covariance function
    is plotted.} 
  \item{dist}{A numeric vector corresponding to the distance at which
    the covariance function should be evaluated. May be missing.}
  \item{xlab,ylab}{The x-axis and y-axis labels. May be missing.}
  \item{\dots}{Several option to be passed to the \code{\link{plot}}
    function.}
}
\value{This function returns the covariance function. Eventually, if
  \code{dist} is given, the covariance function is computed for each
  distance given by \code{dist}. If \code{plot = TRUE}, the covariance
  function is plotted.}
\details{
  Currently, three covariance functions are defined: the Whittle-Matern,
  powered exponential (also known as "stable") and the Cauchy
  models. These covariance functions are defined as follows:

  \describe{
    \item{Whittle-Matern}{\eqn{\rho(h) = sill
    \frac{2^{1-smooth}}{\Gamma(smooth)} \left(\frac{h}{range}
    \right)^{smooth} K_{smooth}\left(frac{h}{range} \right)}{rho(h) =
    sill * 2^(1-smooth) / \Gamma(smooth) * (h/range)^smooth * K_smooth
    (h / range)}}
    \item{Powered Exponential}{\eqn{\rho(h) = sill * \exp \left[-
    \left(\frac{h}{range} \right)^{smooth} \right]}{rho(h) = sill *
    exp[-(h/range)^smooth]}}
    \item{Cauchy}{\eqn{\rho(h) = sill \left[1 + \left(\frac{h}{range}
    \right)^2 \right]^{-smooth}}{rho(h) = sill * [1 +
    (h/range)^2]^(-smooth)}}
  }
  where \eqn{\Gamma}{\Gamma} is the gamma function and
    \eqn{K_{smooth}}{K_smooth} is the modified Bessel function of order
    \eqn{smooth}{smooth}.
}
\examples{
## 1- Calling covariance using fixed covariance parameters
covariance(sill = 1, range = 1, smooth = 0.5, cov.mod = "whitmat")
covariance(sill = 0.5, range = 1, smooth = 0.5, cov.mod = "whitmat",
  dist = seq(0,5, 0.2), plot = FALSE)

## 2- Calling covariance from a fitted model
\dontrun{
require(RandomFields)

##Define the coordinate of each location
n.site <- 30
locations <- matrix(runif(2*n.site, 0, 10), ncol = 2)
colnames(locations) <- c("lon", "lat")

##Simulate a max-stable process - with unit Frechet margins
ms0 <- MaxStableRF(locations[,1], locations[,2], grid=FALSE, model="wh",
                   param=c(0,1,0,3, .5), maxstable="extr",
                   n = 30)
ms1 <- t(ms0)

##Now define the spatial model for the GEV parameters
param.loc <- -10 + 2 * locations[,2]
param.scale <- 5 + 2 * locations[,1] + locations[,2]^2
param.shape <- rep(0.2, n.site)
##Transform the unit Frechet margins to GEV 
for (i in 1:n.site)
  ms1[,i] <- param.scale[i] * (ms1[,i]^param.shape[i] - 1) /
  param.shape[i] + param.loc[i]

##Define a model for the GEV margins to be fitted
##shape ~ 1 stands for a constant GEV shape parameter
##over the region
loc.form <- loc ~ lat
scale.form <- scale ~ lon + I(lat^2)
shape.form <- shape ~ 1

fitted <- fitmaxstab(ms1, locations, "whitmat", loc.form, scale.form,
               shape.form)
covariance(fitted, ylim = c(0, 1))
covariance(sill = 1, range = 3, smooth = .5, cov.mod = "whitmat", add = TRUE, col
               = 3)
title("Whittle-Matern covariance function")
legend("topright", c("Theo.", "Fitted"), lty = 1, col = c(3,1), inset =
           .05)
}
}
\author{Mathieu Ribatet}
\keyword{htest}