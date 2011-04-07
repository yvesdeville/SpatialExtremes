\name{print.copula}
\alias{print.copula}
\title{Printing objects of class copula''}
\description{
  A method for printing object of class ``copula''. 
}
\usage{
\method{print}{copula}(x, digits = max(3, getOption("digits") - 3), ...)
}
\arguments{
  \item{x}{An object of class ``copula. Most often, \code{x} is the
    output of the \code{\link{fitcopula}} function.}
  \item{digits}{The number of digits to be printed.}
  \item{\dots}{Other options to be passed to the \code{\link{print}}
    function.}
}
\value{
  Print several information on screen.
}
\examples{
n.site <- 30
n.obs <- 50

coord <- matrix(runif(2 * n.site, -10, 10), ncol = 2)
colnames(coord) <- c("lon", "lat")

## Generate data from a Gaussian copula model
data <- rcopula(n.obs, coord, "gaussian", "powexp", nugget = 0, range = 4, smooth = 1.2)

## Transform the margins to GEV
locs <- -5 + coord[,"lon"] / 10
scales <- 10 + coord[,"lat"] / 2
shapes <- rep(0.2, n.site)

for (i in 1:n.site)
  data[,i] <- frech2gev(data[,i], locs[i], scales[i], shapes[i])

## Fit a Gaussian copula model

## 1. Define trend surfaces
loc.form <- y ~ lon
scale.form <- y ~ lat
shape.form <- y ~ 1

## 2. Fit
M0 <- fitcopula(data, coord, "gaussian", "powexp", loc.form, scale.form,
                shape.form, nugget = 0)
M0
}
\author{Mathieu Ribatet}
\keyword{print}