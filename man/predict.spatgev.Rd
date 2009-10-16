\name{predict.spatgev}
\alias{predict.spatgev}
\title{Prediction of the GEV parameters}
\description{This function predicts the marginal GEV parameters from a
  fitted ''spatial GEV'' model. 
}
\usage{
\method{predict}{spatgev}(object, newdata, ret.per = NULL, \dots)
}
\arguments{
  \item{object}{An object of class spatgev''. Most often, it will be
    the output of the function \code{\link{fitspatgev}}.}
  \item{newdata}{An optional data frame in which to look for variables
    with which to predict. If omitted, the fitted values are used.}
  \item{ret.per}{Numeric vector giving the return periods for which
    return levels are computed. If \code{NULL} (default), no return
    levels are computed.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\value{'predict.spatgev' produces a vector of predictions or a matrix of
  predictions.
}
\seealso{
  \code{\link{predict}}
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

fitted <- fitspatgev(data, locations, form.loc, form.scale, form.shape)

## 4- GEV parameters estimates at each locations or at ungauged locations
predict(fitted)
ungauged <- data.frame(lon = runif(10, 0, 10), lat = runif(10, 0, 10))
predict(fitted, ungauged)
}
\author{Mathieu Ribatet}
\keyword{htest}