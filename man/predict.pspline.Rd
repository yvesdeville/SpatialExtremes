\name{predict.pspline}
\alias{predict.pspline}
\title{Prediction of smoothing spline with radial basis functions}
\description{This function predicts the response from a fitted
  smoothing spline. 
}
\usage{
\method{predict}{pspline}(object, new.data, \dots)
}
\arguments{
  \item{object}{An object of class ``pspline''. Most often, it will be
    the output of the function \code{\link{rbpspline}}.}
  \item{new.data}{An optional data frame in which to look for variables
    with which to predict. If omitted, the fitted values are used.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\value{'predict.pspline' produces a vector of predictions or a matrix of
  predictions.
}
\seealso{
  \code{\link{predict}}
}
\examples{
## 1- Define a function to approximate
fun <- function(x)
  sin(3 * pi * x)

## 2- Compute the response from fun - and adding a noise
x <- seq(0, 1, length = 200)
y <- fun(x) + rnorm(200, 0, sqrt(0.4))

## 2- Fit a penalized smoothing spline
n.knots <- 30
knots <- quantile(x, prob = 1:n.knots / (n.knots + 2))
fitted <- rbpspline(y, x, knots, degree = 3)

## 3- Prediction from the fitted spline
plot(x, y, pch = 2, col = "red")
plot(fun, from = 0, to = 1, add = TRUE)
pred <- predict(fitted)
lines(pred[,1], pred[,2], col = "blue", pch = 3)
}
\author{Mathieu Ribatet}
\keyword{htest}