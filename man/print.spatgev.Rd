\name{print}
\alias{print.spatgev}
\title{Printing objects of class ''spatgev''}
\description{
  A method for printing object of class ''spatgev''. 
}
\usage{
\method{print}{spatgev}(x, digits = max(3, getOption("digits") - 3),
\dots)
}
\arguments{
  \item{x}{An object of class ''spatgev''. Most often, \code{x} is the
    output of the \code{\link{fitspatgev}} function.}
  \item{digits}{The number of digits to be printed.}
  \item{\dots}{Other options to be passed to the \code{\link{print}}
    function.}
}
\value{
  Print several information on screen.
}
\examples{
## 1- Simulate a max-stable random field
require(RandomFields)
n.site <- 35
locations <- matrix(runif(2*n.site, 0, 10), ncol = 2)
colnames(locations) <- c("lon", "lat")
ms0 <- MaxStableRF(locations[,1], locations[,2], grid=FALSE, model="wh",
                   param=c(0,1,0,3, .5), maxstable="extr",
                   n = 50)
## 2- Transformation to non unit Frechet margins
ms1 <- t(ms0)
param.loc <- -10 + 2 * locations[,2]
param.scale <- 5 + 2 * locations[,1]
param.shape <- rep(0.2, n.site)
for (i in 1:n.site)
  ms1[,i] <- param.scale[i] * (ms1[,i]^param.shape[i] - 1) /
  param.shape[i] + param.loc[i]

## 3- Fit a ''spatial GEV'' mdoel to data with the following models for
##    the GEV parameters
form.loc <- loc ~ lat
form.scale <- scale ~ lon
form.shape <- shape ~ 1

fitspatgev(ms1, locations, form.loc, form.scale, form.shape)
}
\author{Mathieu Ribatet}
\keyword{print}