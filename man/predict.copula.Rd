\name{predict.copula}
\alias{predict.copula}
\title{Prediction of the copula marginal parameters}
\description{This function predicts the marginal GEV parameters from a
  fitted copula based model. 
}
\usage{
\method{predict}{copula}(object, newdata, ret.per = NULL, std.err =
TRUE, \dots)
}
\arguments{
  \item{object}{An object of class ``copula''. Most often, it will be
    the output of the function \code{\link{fitcopula}}.}
  \item{newdata}{An optional data frame in which to look for variables
    with which to predict. If omitted, the fitted values are used.}
  \item{ret.per}{Numeric vector giving the return periods for which
    return levels are computed. If \code{NULL} (default), no return
    levels are computed.}
  \item{std.err}{Logical. If \code{TRUE} (default), standard erros will
    be computed if possible.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\value{'predict.copula' produces a vector of predictions or a matrix of
  predictions.
}
\seealso{
  \code{\link{predict}}
}
\examples{
## 1- Simulate a max-stable random field
n.site <- 10
locations <- matrix(runif(2*n.site, 0, 10), ncol = 2)
colnames(locations) <- c("lon", "lat")

data <- rmaxstab(25, locations, cov.mod = "whitmat", nugget = 0, range = 30,
smooth = 0.5)

## 2- Transformation to non unit Frechet margins
param.loc <- -10 + 2 * locations[,2]
param.scale <- 5 + 2 * locations[,1]
param.shape <- rep(0.2, n.site)

for (i in 1:n.site)
  data[,i] <- frech2gev(data[,i], param.loc[i], param.scale[i],
param.shape[i])

## 3- Fit a copula based model with the following model for
##    the GEV parameters
form.loc <- loc ~ lat
form.scale <- scale ~ lon
form.shape <- shape ~ 1

M0 <- fitcopula(data, locations, "student", "whitmat", loc.form = form.loc,
                        scale.form = form.scale, shape.form =
                        form.shape, nugget = 0)

## 4- GEV parameters estimates at each locations or at ungauged locations
predict(M0)
ungauged <- data.frame(lon = runif(10, 0, 10), lat = runif(10, 0, 10))
predict(M0, ungauged)
}
\author{Mathieu Ribatet}
\keyword{htest}