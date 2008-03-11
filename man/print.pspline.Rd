\name{print}
\alias{print.pspline}
\title{Printing objects of class ``pspline''}
\description{
  A method for printing object of class ``pspline''. 
}
\usage{
\method{print}{pspline}(x, \dots)
}
\arguments{
  \item{x}{An object of class ``pspline''. Most often, \code{x} is the
    output of the \code{\link{rbpspline}} function.}
  \item{\dots}{Other options to be passed to the \code{\link{print}}
    function.}
}
\value{
  Print several information on screen.
}
\examples{
n <- 200
x <- runif(n)
fun <- function(x) sin(3 * pi * x)
y <- fun(x) + rnorm(n, 0, sqrt(0.4))
knots <- quantile(x, prob = 1:31 / 32)
fitted <- rbpspline(y, x, knots = knots, degree = 3)
fitted
}
\author{Mathieu Ribatet}
\keyword{print}