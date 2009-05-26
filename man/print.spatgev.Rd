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
n.site <- 35
locations <- matrix(runif(2*n.site, 0, 10), ncol = 2)
colnames(locations) <- c("lon", "lat")

data <- rmaxstab(50, locations, cov.mod = "whitmat", sill = 1, range = 3,
smooth = 0.5)

## 2- Transformation to non unit Frechet margins
param.loc <- -10 + 2 * locations[,2]
param.scale <- 5 + 2 * locations[,1]
param.shape <- rep(0.2, n.site)

for (i in 1:n.site)
  data[,i] <- frech2gev(data[,i], param.loc[i], param.scale[i],
param.shape[i])

## 3- Fit a ''spatial GEV'' mdoel to data with the following models for
##    the GEV parameters
form.loc <- loc ~ lat
form.scale <- scale ~ lon
form.shape <- shape ~ 1

fitspatgev(data, locations, form.loc, form.scale, form.shape)
}
\author{Mathieu Ribatet}
\keyword{print}