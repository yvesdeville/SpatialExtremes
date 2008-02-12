\name{fitmaxstab}
\alias{fitmaxstab}
\title{Fits a max-stable process to data}
\description{This function fits max-stable processes to data using
  pairwise likelihood. Two max-stable characterisations are available:
  the Smith and Schlather representations.
}
\usage{fitmaxstab(data,  coord, cov.mod, loc.form, scale.form,
shape.form, marg.cov = NULL, iso = FALSE, \dots, fit.marge = FALSE, warn
= TRUE, method = "Nelder", start, control = list(), std.err.type =
"score", corr = FALSE)}
\arguments{
  \item{data}{A matrix representing the data. Each column corresponds to
    one location.}
  \item{coord}{A matrix that gives the coordinates of each
    location. Each row corresponds to one location.}
  \item{cov.mod}{A character string corresponding the the covariance
    model in the max-stable representation. Must be one of "gauss" for
    the Smith's model; or "whitmat", "cauchy" or "powexp" for
    the Whittle-Matern, the Cauchy and the Powered Exponential
    covariance family with the Schlather's model. These 3 last
    covariance models could have an additional prefix "i" or "g" for the
    mixture of independent Frechet field and Schlather's model or the
    geometric gaussian model.} 
  \item{loc.form, scale.form, shape.form}{R formulas defining the
    spatial linear model for the GEV parameters. May be missing. See
    section details.}
  \item{marg.cov}{Matrix with named columns giving additional covariates
    for the GEV parameters. If \code{NULL}, no extra covariates are
    used.}
  \item{iso}{Logical. If \code{TRUE} an isotropic model is fitted to
    data. Otherwise (default), anisotropy is allowed. Currently, this is
    only implemented for the Smith's model.}
  \item{\dots}{Several arguments to be passed to the \code{\link{optim}}
    function. See section details.}
  \item{fit.marge}{Logical. If \code{TRUE}, the GEV parameters
    are estimated ``point-wisely'' or using the formulas given by
    \code{loc.form}, \code{scale.form} and \code{shape.form}. If
    \code{FALSE}, observations are supposed to be unit Frechet
    distributed. Note that when formulas are given, \code{fit.marge} is
    automatically set to \code{TRUE}.}
  \item{warn}{Logical. If \code{TRUE} (default), users are warned if
    the log-likelihood is infinite at starting values and/or problems
    arised while computing the standard errors.}
  \item{method}{The method used for the numerical optimisation
    procedure. Must be one of \code{BFGS}, \code{Nelder-Mead},
    \code{CG}, \code{L-BFGS-B}, \code{SANN} or \code{nlm}. See \code{\link{optim}}
    for details. Please note that passing \code{nlm} will use the
    \code{\link{nlm}} function instead of \code{\link{optim}}.}
  \item{start}{A named list giving the initial values for the
    parameters over which the pairwise likelihood is to be minimized. If
    \code{start} is omitted the routine attempts to find good starting
    values - but might fail.}
  \item{control}{A list giving the control parameters to be passed to
    the \code{\link{optim}} function.}
  \item{std.err.type}{Character string. Must be one of "none", "score",
    "grad". If "none", standard errors are not computed. Otherwise,
    standard errors are estimated using the sandwich estimates - see
    section Details.}
  \item{corr}{Logical. If \code{TRUE} (non default), the asymptotic
    correlation matrix is computed.}
}
\value{
  This function returns a object of class \code{maxstab}. Such objects
  are list with components:
  \item{fitted.values}{A vector containing the estimated parameters.}
  \item{std.err}{A vector containing the standard errors.}
  \item{fixed}{A vector containing the parameters of the model that
    have been held fixed.}
  \item{param}{A vector containing all parameters (optimised and fixed).}
  \item{deviance}{The (pairwise) deviance at the maximum pairwise
    likelihood estimates.}
  \item{corr}{The correlation matrix.}
  \item{convergence, counts, message}{Components taken from the
    list returned by \code{\link{optim}} - for the \code{mle} method.}
  \item{data}{The data analysed.}
  \item{model}{The max-stable characterisation used.}
  \item{fit.marge}{A logical that specifies if the GEV margins were
    estimated.}
  \item{cov.fun}{The estimated covariance function - for the Schlather
    model only.}
  \item{extCoeff}{The estimated extremal coefficient function.}
  \item{cov.mod}{The covariance model for the spatial structure.}
}
\details{
  The function is a wrapper to more specialised functions
  (\code{schlatherfull, schlatherlm, smithfull, smithlm}).

  As spatial data often deal with a large number of locations, it is
  impossible to write analytically the joint distribution. Consequently,
  the fitting procedure substitutes the "full likelihood" for the
  pairwise likelihood.

  Let define \eqn{L_{i,j}(x_{i,j}, \theta)}{L_{i,j}(x_{i,j}, theta)} the
  likelihood for site \eqn{i}{i} and \eqn{j}{j}, where \eqn{i = 1,
  \dots, N-1}{i = 1, \dots, N-1}, \eqn{j = i+1, \dots, N}{j=i+1, \dots, 
  N}, \eqn{N}{N} is the number of site within the region and
  \eqn{x_{i,j}}{x_{i,j}} are the joint observations for site \eqn{i}{i}
  and \eqn{j}{j}. Then the pairwise likelihood
  \eqn{PL(\theta)}{PL(theta)} is defined by:

  \deqn{\ell_P = \log PL(\theta) = \sum_{i = 1}^{N-1} \sum_{j=i+1}^{N} \log
  L_{i,j} (x_{i,j}, \theta)}{llik_P = log PL(theta) = sum_{i = 1}^{N-1}
  sum_{j=i+1}^{N} log L_{i,j} (x_{i,j}, theta)}

  As pairwise likelihood is an approximation of the ``full likelihood'',
  standard errors cannot be computed directly by the inverse of the
  Fisher information matrix. Instead, a sandwich estimate must be used
  to account for model mispecification e.g.

  \deqn{\hat{\theta} \sim N(\theta, H^{-1} J H^{-1})}{hat(\theta) ~
    N(\theta, H^{-1} J H^{-1})}
  where \eqn{H} is the Fisher information matrix (computed from the
    misspecified model) and \eqn{J} is the variance of the score
    function.

  \eqn{H} is easily estimated using the observed Hessian matrix given by
  the \code{\link{optim}} function. Estimation of \eqn{J} is much more
  difficult. Currently, we propose two different strategies to estimate
  \eqn{J}.

  \describe{
    \item{\code{grad}}{\eqn{J} is estimated from the gradient e.g. \eqn{J
	= \nabla \ell_P \nabla \ell_P^T}{J = grad(llik_P) grad(llik_P)^T}}
    \item{\code{score}}{\eqn{J} is estimated directly from the variance of
      the observed score function.}
  }
}
\references{
  Cox, D. R. and Reid, N. (2004) A note on pseudo-likelihood constructed
  from marginal densities. \emph{Biometrika} \bold{91}: 729--737.

  Schlather, M. (2002) Models for Stationary Max-Stable Random
  Fields. \emph{Extremes} \bold{5}:1,  33--44.

  Smith, R. L. (1990) Max-stable processes and spatial
  extremes. Unpublished.
}
\examples{
require(RandomFields)

##Define the coordinate of each location
n.site <- 30
locations <- matrix(runif(2*n.site, 0, 10), ncol = 2)
colnames(locations) <- c("lon", "lat")

##Simulate a max-stable process - with unit Frechet margins
ms0 <- MaxStableRF(locations[,1], locations[,2], grid=FALSE, model="wh",
                   param=c(0,1,0,30, .5), maxstable="extr",
                   n = 40)
ms0 <- t(ms0)
ms1 <- ms0

##Now define the spatial model for the GEV parameters
param.loc <- -10 + 2 * locations[,2]
param.scale <- 5 + 2 * locations[,1] + locations[,2]^2
param.shape <- rep(0.2, n.site)

##Transform the unit Frechet margins to GEV 
for (i in 1:n.site)
  ms1[,i] <- param.scale[i] * (ms1[,i]^param.shape[i] - 1) /
  param.shape[i] + param.loc[i]

##Define a model for the GEV margins to be fitted
##shape ~ 1 stands for the GEV shape parameter is constant
##over the region
loc.form <- loc ~ lat
scale.form <- scale ~ lon + I(lat^2)
shape.form <- shape ~ 1

##Fit a max-stable process using the Schlather's model
fitmaxstab(ms1, locations, "whitmat", loc.form, scale.form,
           shape.form)

## Model without any spatial structure for the GEV parameters
## Be careful this could be *REALLY* time consuming
\dontrun{fitmaxstab(ms1, locations, "whitmat")}

##  Fixing the smooth parameter of the Whittle-Matern family
##  to 0.5 - e.g. considering exponential family. We suppose the date
##  are unit Frechet here.
\dontrun{fitmaxstab(ms0, locations, "whitmat", smooth = 0.5, fit.marge = FALSE)}

##  Fitting a penalized smoothing splines for the margins with the
##     Smith's model
\dontrun{
sigma <- matrix(c(100, 25, 25, 220),ncol = 2)
sigma.inv <- solve(sigma)
sqrtCinv <- t(chol(sigma.inv))
model <- list(list(model = "gauss", var = 1, aniso = sqrtCinv / 2))
ms0 <- MaxStableRF(locations[,1], locations[,2], grid=FALSE, model=model,
                   maxstable="Bool", n = 40)
ms1 <- t(ms0)
}
##     And transform it to ordinary GEV margins with a non-linear
##     function
\dontrun{
fun <- function(x)
  2 * sin(pi * x / 4) + 10
fun2 <- function(x)
  (fun(x) - 7 ) / 15

param.loc <- fun(locations[,2])
param.scale <- fun(locations[,2])
param.shape <- fun2(locations[,1])
}

##Transformation from unit Frechet to common GEV margins
\dontrun{
for (i in 1:n.site)
  ms1[,i] <- param.scale[i] * (ms1[,i]^param.shape[i] - 1) / param.shape[i] + param.loc[i]
}

##Defining the knots, penalty, degree for the splines
\dontrun{
n.knots <- 5
knots <- quantile(locations[,2], prob = 1:n.knots/(n.knots+1))
knots2 <- quantile(locations[,1], prob = 1:n.knots/(n.knots+1))
}

##Be careful the choice of the penalty (i.e. the smoothing parameter)
##may strongly affect the result Here we use p-splines for each GEV
##parameter - so it's really CPU demanding but one can use 1 p-spline
##and 2 linear models.
##A simple linear model will be clearly faster...
\dontrun{
loc.form <- y ~ rb(lat, knots = knots, degree = 3, penalty = .5)
scale.form <- y ~ rb(lat, knots = knots, degree = 3, penalty = .5)
shape.form <- y ~ rb(lon, knots = knots2, degree = 3, penalty = .5)

fitted <- fitmaxstab(ms1, locations, "gauss", loc.form, scale.form, shape.form,
                     control = list(ndeps = rep(1e-6, 24), trace = 10),
                     std.err.type = "none", method = "BFGS")
fitted
op <- par(mfrow=c(1,3))
plot(locations[,2], param.loc, col = 2, ylim = c(7, 14),
     ylab = "location parameter", xlab = "latitude")
plot(fun, from = 0, to = 10, add = TRUE, col = 2)
points(locations[,2], predict(fitted)[,"loc"], col = "blue", pch = 5)
new.data <- cbind(lon = seq(0, 10, length = 100), lat = seq(0, 10, length = 100))
lines(new.data[,1], predict(fitted, new.data)[,"loc"], col = "blue")
legend("topleft", c("true values", "predict. values", "true curve", "predict. curve"),
       col = c("red", "blue", "red", "blue"), pch = c(1, 5, NA, NA), inset = 0.05,
       lty = c(0, 0, 1, 1), ncol = 2)

plot(locations[,2], param.scale, col = 2, ylim = c(7, 14),
     ylab = "scale parameter", xlab = "latitude")
plot(fun, from = 0, to = 10, add = TRUE, col = 2)
points(locations[,2], predict(fitted)[,"scale"], col = "blue", pch = 5)
lines(new.data[,1], predict(fitted, new.data)[,"scale"], col = "blue")
legend("topleft", c("true values", "predict. values", "true curve", "predict. curve"),
       col = c("red", "blue", "red", "blue"), pch = c(1, 5, NA, NA), inset = 0.05,
       lty = c(0, 0, 1, 1), ncol = 2)

plot(locations[,1], param.shape, col = 2,
     ylab = "shape parameter", xlab = "longitude")
plot(fun2, from = 0, to = 10, add = TRUE, col = 2)
points(locations[,1], predict(fitted)[,"shape"], col = "blue", pch = 5)
lines(new.data[,1], predict(fitted, new.data)[,"shape"], col = "blue")
legend("topleft", c("true values", "predict. values", "true curve", "predict. curve"),
       col = c("red", "blue", "red", "blue"), pch = c(1, 5, NA, NA), inset = 0.05,
       lty = c(0, 0, 1, 1), ncol = 2)
par(op)
}
}
\author{Mathieu Ribatet}
\keyword{htest}
