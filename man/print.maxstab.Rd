\name{print}
\alias{print.maxstab}
\title{Printing objects of class ``maxstab''}
\description{
  A method for printing object of class ``maxstab''. 
}
\usage{
\method{print}{maxstab}(x, digits = max(3, getOption("digits") - 3),
\dots)
}
\arguments{
  \item{x}{An object of class ``maxstab''. Most often, \code{x} is the
    output of the \code{\link{fitmaxstab}} function.}
  \item{digits}{The number of digits to be printed.}
  \item{\dots}{Other options to be passed to the \code{\link{print}}
    function.}
}
\value{
  Print several information on screen.
}
\examples{
require(RandomFields)

##Define the coordinates of each location
n.site <- 30
locations <- matrix(5 + rnorm(2*n.site, sd = sqrt(2)), ncol = 2)

##Simulate a max-stable process - with unit Frechet margins
ms0 <- MaxStableRF(locations[,1], locations[,2], grid=FALSE, model="wh",
                   param=c(0,1,0,3, .5), maxstable="extr",
                   n = 30)
ms1 <- t(ms0)
\dontrun{
fit <- fitmaxstab(ms1, locations, "schlather")
fit
}
}
\author{Mathieu Ribatet}
\keyword{print}