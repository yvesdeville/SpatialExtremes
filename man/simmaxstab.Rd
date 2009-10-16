\name{rmaxstab}
\alias{rmaxstab}

\title{Simulation of Max-Stable Random Fields}
\description{This functions generates realisation from a max-stable
  random field.}

\usage{
rmaxstab(n, coord, cov.mod = "gauss", grid = FALSE, control =
list(), \dots)
}

\arguments{
  \item{n}{Integer. The number of observations.}
  \item{coord}{A vector or matrix that gives the coordinates of each
    location. Each row corresponds to one location - if any.}
  \item{cov.mod}{A character string that gives the max-stable
    model. This must be one of "gauss" for the Smith model or "whitmat",
    "cauchy", "powexp" and "bessel" for the Schlather model with the given
    correlation family.}
  \item{grid}{Logical. Does the coordinates represent grid points?}
  \item{control}{A named list with arguments 'nlines' (number of lines
    of the TBM simulation), 'method' the name of the simulation
    method - must be one of 'exact' or 'tbm', and 'uBound' the uniform
  upper bound. See details.}
  \item{\dots}{The parameters of the max-stable model. See details.}
}
\details{
  Users must supply the parameters for the max-stable model. For the
  Schlather model, users should supply the "sill", "range" and "smooth"
  parameter values. For the Smith model, if \code{coord} is univariate
  you must specify \code{var}, otherwise users should supply the
  covariance parameters i.e. parameters with names such as \code{cov11},
  \code{cov12}, \ldots

  Here are the details for each allowed components of 'control'. If
  'method' is \code{NULL} (default), the function tries to find the 
  most appropriate simulation technique. If 'nlines' is \code{NULL}, it
  is set to 1000. If 'uBound' is \code{NULL}, it is set to reasonable
  values - for example 3.5 for the Schalther model.
}
\value{
  A matrix containing observations from the required max-stable
  model. Each column represents one stations. If \code{grid = TRUE}, the
  function returns an array of dimension nrow(coord) x nrow(coord) x n.
}
\references{
  Schlather, M. (2002) Models for Stationary Max-Stable Random
  Fields. \emph{Extremes} \bold{5}:1,  33--44.
  
  Smith, R. L. (1990) Max-stable processes and spatial
  extremes. Unpublished manuscript.
}
\author{Mathieu Ribatet}
\seealso{\code{fitmaxstab}}
\examples{
## 1. Smith's model
set.seed(8)
x <- seq(0, 10, length = 200)
coord <- cbind(x, x)
data <- rmaxstab(1, coord, cov11 = 9/8, cov12 = 0, cov22 = 9/8, grid = TRUE)
##We change the margins for visibility
filled.contour(x, x, sqrt(data[,,1]))

## 2. Schlather's model
coord <- matrix(runif(100, 0, 15), ncol = 2)
data <- rmaxstab(100, coord, cov.mod = "whitmat", sill = 1, range = 10,
smooth = 1)
}
\keyword{distribution}
