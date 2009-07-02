\name{fitspatgev}
\alias{fitspatgev}
\title{MLE for a spatial GEV model}
\description{
  This function derives the MLE of a spatial GEV model.
}
\usage{
fitspatgev(data, covariables, loc.form, scale.form, shape.form, ...,
start, control = list(maxit = 10000), method = "Nelder", std.err.type =
"score", warn = TRUE)
}
\arguments{
  \item{data}{A matrix representing the data. Each column corresponds to
    one location.}
  \item{covariables}{Matrix with named columns giving required
    covariates for the GEV parameter models.}
  \item{loc.form, scale.form, shape.form}{R formulas defining the
    spatial  models for the GEV parameters. See section details.}
  \item{start}{A named list giving the initial values for the
    parameters over which the pairwise likelihood is to be minimized. If
    \code{start} is omitted the routine attempts to find good starting
    values - but might fail.}
  \item{\dots}{Several arguments to be passed to the
    \code{\link{optim}} functions. See section details.}
  \item{control}{The control argument to be passed to the
    \code{\link{optim}} function.}
  \item{method}{The method used for the numerical optimisation
    procedure. Must be one of \code{BFGS}, \code{Nelder-Mead},
    \code{CG}, \code{L-BFGS-B} or \code{SANN}. See \code{\link{optim}}
    for details.}
  \item{std.err.type}{Character string. Must be one of "score", "grad"
    or "none". If none, no standard errors are computed.}
  \item{warn}{Logical. If \code{TRUE} (default), users will be warned if
  the starting values lie in a zero density region.}
}
\details{
  A kind of "spatial" GEV model can be defined by using response
  surfaces for the GEV parameters. For instance, the GEV location
  parameters are defined through the following equation:

  \deqn{\mu = X_\mu \beta_\mu}{\mu = X_\mu \beta_\mu}
  where \eqn{X_\mu}{X_\mu} is the design matrix and
  \eqn{\beta_\mu}{\beta_\mu} is the vector parameter to be
  estimated. The GEV scale and shape parameters are defined accordingly
  to the above equation.

  The log-likelihood for the GEV spatial model is consequently defined
  as follows:

  \deqn{\ell(\beta) = \sum_{i=1}^{n.site} \sum_{j=1}^{n.obs} \log
    f(y_{i,j}; \theta_i)}{llik(\beta) = sum_(i=1)^(n.site)
    \sum_(j=1)^(n.obs) log f(y_(i,j);theta_i)}
  where \eqn{\theta_i}{theta_i} is the vector of the GEV parameters for
  the \eqn{i}-th site.

  Most often, there will be some dependence between stations. However,
  it can be seen from the log-likelihood definition that we supposed
  that the stations are mutually independent. Consequently, to get
  reliable standard error estimates, these standard errors are estimated
  with their sandwich estimates.  
}
\value{
  An object of class \code{spatgev}. Namely, this is a list with the
  following arguments:
  \item{fitted.values}{The parameter estimates.}
  \item{param}{All the parameters e.g. parameter estimates and fixed
    parameters.}
  \item{std.err}{The standard errors.}
  \item{var.cov}{The asymptotic MLE variance covariance matrix.}
  \item{counts,message,convergence}{Some information about the
    optimization procedure.}
  \item{logLik,deviance}{The log-likelihood and deviance values.}
  \item{loc.form, scale.form, shape.form}{The formulas defining the
    spatial models for the GEV parameters.}
  \item{covariables}{The covariables used for the spatial models.}
  \item{ihessian}{The inverse of the Hessian matrix of the negative
    log-likelihood.}
  \item{jacobian}{The variance covariance matrix of the score.}
}
\author{Mathieu Ribatet}
\examples{
## 1- Simulate a max-stable random field
require(RandomFields)
n.site <- 35
locations <- matrix(runif(2*n.site, 0, 10), ncol = 2)
colnames(locations) <- c("lon", "lat")
ms0 <- MaxStableRF(locations[,1], locations[,2], grid=FALSE, model="wh",
                   param=c(0,1,0,3, .5), maxstable="extr",
                   n = 50)
## 2- Transformation to non unit Frechet margins
ms1 <- t(ms0)
param.loc <- -10 + 2 * locations[,2]
param.scale <- 5 + 2 * locations[,1]
param.shape <- rep(0.2, n.site)
for (i in 1:n.site)
  ms1[,i] <- param.scale[i] * (ms1[,i]^param.shape[i] - 1) /
  param.shape[i] + param.loc[i]

## 3- Fit a ''spatial GEV'' mdoel to data with the following models for
##    the GEV parameters
form.loc <- loc ~ lat
form.scale <- scale ~ lon
form.shape <- shape ~ 1

fitspatgev(ms1, locations, form.loc, form.scale, form.shape)
}
\keyword{htest}

