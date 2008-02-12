\name{fitcovmat}
\alias{fitcovmat}
\title{Estimates the covariance matrix for the Smith's model}
\description{
  Estimates the covariance matrix for the Smith's model using
  non-parametric estimates of the pairwise extremal coefficients.
}
\usage{
fitcovmat(data, coord, marge = "mle", iso = FALSE, start, \dots)
}
\arguments{
  \item{data}{A matrix representing the data. Each column corresponds to
    one location.}
  \item{coord}{A matrix that gives the coordinates of each
    location. Each row corresponds to one location.}
   \item{marge}{Character string specifying how margins are transformed
    to unit Frechet. Must be one of "emp", "frech" or "mle" - see
    function \code{\link{fitextcoeff}}.}
  \item{iso}{Logical. If \code{TRUE}, isotropy is supposed. Otherwise
    (default), anisotropy is allowed.}
  \item{start}{A named list giving the initial values for the
    parameters over which the weighted sum of square is to be
    minimized. If \code{start} is omitted the routine attempts to find
    good starting values.}
  \item{\dots}{Optional arguments to be passed to the \code{\link{nlm}}
    function.}
}
\details{
  The fitting procedure is based on weighted least squares. More
  precisely, the fitting criteria is to minimize:
  \deqn{\sum_{i,j} \left(\frac{\tilde{\theta}_{i,j} -
      \hat{\theta}_{i,j}}{s_{i,j}}\right)^2}{\sum_{i,j}
    [(\theta_{i,j}^+ - \theta_{i,j}^*) / s_{i,j}]^2}
  where \eqn{\tilde{\theta}_{i,j}}{\theta_{i,j}^+} is a non
  parametric estimate of the extremal coefficient related to location
  \code{i} and \code{j}, \eqn{\hat{\theta}_{i,j}}{\theta_{i,j}^*} is
  the fitted extremal coefficient derived from the Smith's model and
  \eqn{s_{i,j}}{s_{i,j}} are the standard errors related to the
  estimates \eqn{\tilde{\theta}_{i,j}}{\theta_{i,j}^+}.
}
\value{
  An object of class maxstab.
}
\references{
  Smith, R. L. (1990) Max-stable processes and spatial
  extremes. Unpublished manuscript.
  }
\author{Mathieu Ribatet}
\seealso{\code{\link{fitcovariance}}, \code{\link{fitmaxstab}},
  \code{\link{fitextcoeff}}}
\examples{
require(RandomFields)
n.site <- 50
n.obs <- 100

locations <- matrix(runif(2*n.site, 0, 40), ncol = 2)
colnames(locations) <- c("lon", "lat")

sigma <- matrix(c(200, 0, 0, 200),ncol = 2)
sigma.inv <- solve(sigma)
sqrtCinv <- t(chol(sigma.inv))

model <- list(list(model = "gauss", var = 1, aniso = sqrtCinv / 2))

## Simulate a max-stable process - with unit Frechet margins
ms1 <- MaxStableRF(locations, model = model, maxstable = "Bool", grid = FALSE, n = n.obs)
ms1 <- t(ms1)

fitcovmat(ms1, locations, marge = "emp")

##Force an isotropic model
fitcovmat(ms1, locations, marge = "emp", iso = TRUE)
}
\keyword{htest}

