\name{TIC}
\alias{TIC}
\alias{TIC.maxstab}
\alias{TIC.spatgev}
\title{Takeuchi's information criterion}
\description{
  Computes a "generalization" of the Takeuchi's information criterion
  when the model is miss-specified.
}
\usage{
\method{TIC}{maxstab}(object, \dots)
\method{TIC}{spatgev}(object, \dots)
}
\arguments{
  \item{object}{An object of class \code{maxstab} or
    \code{spatgev}. Often, it will be the output of the
    \code{\link{fitmaxstab}} or \code{\link{fitspatgev}} function.}
  \item{\dots}{Additional objects of class \code{maxstab} or
    \code{spatgev} for which TIC should be computed.}
}
\value{
  Numeric.
}
\details{
  TIC is like AIC so that when comparing models one wants to get the
  lowest TIC score.
}
\references{
  Varin, C. and Vidoni, P. (2005) A note on composite likelihood
  inference and model selection. Biometrika \bold{92}(3):519--528.
}
\seealso{
  \code{\link{fitmaxstab}}, \code{\link{AIC}}
}
\examples{
require(RandomFields)

##Define the coordinate of each location
n.site <- 50
locations <- matrix(runif(2*n.site, 0, 100), ncol = 2)
colnames(locations) <- c("lon", "lat")

##Simulate a max-stable process - with unit Frechet margins
ms0 <- MaxStableRF(locations[,1], locations[,2], grid=FALSE, model="stable",
                   param=c(0,1,.2,30, .5), maxstable="extr",
                   n = 40)
ms0 <- t(ms0)
M0 <- fitmaxstab(ms0, locations, "powexp", std.err.type = "score",
                     fit.marge = FALSE)
M1 <- fitmaxstab(ms0, locations, "cauchy", std.err.type = "score",
                     fit.marge = FALSE)
TIC(M0, M1)
}
\author{Mathieu Ribatet}
\keyword{htest}
