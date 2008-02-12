\name{qqextcoeff}
\alias{qqextcoeff}
\title{QQ-plot for the extremal coefficient}
\description{
  This function compares the extremal coefficients estimated from a
  fitted max-stable process to the ones estimated semi-parametrically.
}
\usage{
qqextcoeff(fitted, estim = "ST", marge = "emp", xlab = "Semi-Empirical",
ylab = "Model", ...)
}
\arguments{
  \item{fitted}{An object of class \code{maxstab}. Most often, this will
  be the output of \code{\link{fitmaxstab}}, \code{\link{fitcovmat}} or
  \code{\link{fitcovariance}}.}
  \item{estim,marge}{The \code{estim} and \code{marge} options to be
    passed to the \code{\link{fitextcoeff}} function.}
  \item{xlab,ylab}{The x and y-axis labels.}
  \item{\dots}{Optional arguments to be passed to the \code{\link{plot}}
  function.}
}
\value{
  A QQ-plot.
}
\references{
  Schlather, M. (2002) Models for Stationary Max-Stable Random
  Fields. \emph{Extremes} \bold{5}:1,  33--44.
  
  Schlather, M. and Tawn, J. A. (2003) A dependence measure for
  multivariate and spatial extreme values: Properties and
  inference. Biometrika \bold{90}(1):139--156.

  Smith, R. L. (1990) Max-stable processes and spatial
  extremes. Unpublished manuscript.
}
\author{Mathieu Ribatet}
\seealso{\code{\link{fitmaxstab}}, \code{\link{fitextcoeff}}}
\examples{
require(RandomFields)

##Define the coordinate of each location
n.site <- 30
locations <- matrix(runif(2*n.site, 0, 10), ncol = 2)
colnames(locations) <- c("lon", "lat")

##Simulate a max-stable process - with unit Frechet margins
sigma <- matrix(c(100, 25, 25, 220),ncol = 2)
sigma.inv <- solve(sigma)
sqrtCinv <- t(chol(sigma.inv))
model <- list(list(model = "gauss", var = 1, aniso = sqrtCinv / 2))
ms0 <- MaxStableRF(locations[,1], locations[,2], grid=FALSE, model=model,
                   maxstable="Bool", n = 40)
ms0 <- t(ms0)

fitted <- fitmaxstab(ms0, locations, "gauss", fit.marge = FALSE)

qqextcoeff(fitted)
}
\keyword{hplot}
