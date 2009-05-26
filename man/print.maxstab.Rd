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
##Define the coordinates of each location
n.site <- 30
locations <- matrix(5 + rnorm(2*n.site, sd = sqrt(2)), ncol = 2)

##Simulate a max-stable process - with unit Frechet margins
data <- rmaxstab(30, locations, cov.mod = "whitmat", sill = 1, range = 3,
smooth = 0.5)

fit <- fitmaxstab(data, locations, "whitmat")
fit
}
\author{Mathieu Ribatet}
\keyword{print}