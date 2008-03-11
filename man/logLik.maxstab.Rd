\name{logLik}
\alias{logLik}
\alias{logLik.maxstab}
\title{Extracts Log-Likelihood}
\description{Extract the pairwise log-likelihood for objects of class
  ``maxstab''}
\usage{
\method{logLik}{maxstab}(object, \dots)
}
\arguments{
  \item{object}{An object of class ``maxstab''. Most often this will be
    the output of the \code{\link{fitmaxstab}} function.}
  \item{\dots}{Other arguments to be passed to the \code{\link{logLik}}
    function.}
}
\value{
  Standard \code{logLik} object (see \code{\link{logLik}}) except that
  this is the pairwise log-likelihood! 
}
\seealso{
  \code{\link{logLik}}
}
\examples{
require(RandomFields)

##Define the coordinates of each location
n.site <- 30
locations <- matrix(5 + runif(2*n.site, 0, 10), ncol = 2)

##Simulate a max-stable process - with unit Frechet margins
ms0 <- MaxStableRF(locations[,1], locations[,2], grid=FALSE, model="wh",
                   param=c(0,1,0,3, .5), maxstable="extr",
                   n = 30)
ms1 <- t(ms0)
\dontrun{
fit <- fitmaxstab(ms1, locations, "schlather")
logLik(fit)
}
}
\author{Mathieu Ribatet}
\keyword{models}