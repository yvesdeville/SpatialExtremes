\name{covariance}
\alias{covariance}
\title{Defines and computes covariance functions}
\description{This function defines and computes several covariance
  function either from a fitted ``max-stable'' model; either by
  specifying directly the covariance parameters.
}
\usage{
covariance(fitted, sill, range, smooth, cov.mod = "whitmat", plot = TRUE,
dist, xlab, ylab, ...)
}
\arguments{
  \item{fitted}{An object of class ``maxstab''. Most often this will be
    the output of the \code{\link{fitmaxstab}} function. May be missing
    if \code{sill}, \code{range}, \code{smooth} and \code{cov.mod} are
    given.}
  \item{sill,range,smooth}{The sill, scale and smooth parameters for
    the covariance function. May be missing if \code{fitted} is given.}
  \item{cov.mod}{Character string. The name of the covariance
    model. Must be one of "whitmat", "cauchy", "powexp" or "bessel" for
    the Whittle-Matern, Cauchy, Powered Exponential and Bessel
    models. May be missing if \code{fitted} is given.}
  \item{plot}{Logical. If \code{TRUE} (default) the covariance function
    is plotted.} 
  \item{dist}{A numeric vector corresponding to the distance at which
    the covariance function should be evaluated. May be missing.}
  \item{xlab,ylab}{The x-axis and y-axis labels. May be missing.}
  \item{\dots}{Several option to be passed to the \code{\link{plot}}
    function.}
}
\value{This function returns the covariance function. Eventually, if
  \code{dist} is given, the covariance function is computed for each
  distance given by \code{dist}. If \code{plot = TRUE}, the covariance
  function is plotted.}
\details{
  Currently, four covariance functions are defined: the Whittle-Matern,
  powered exponential (also known as "stable"), Cauchy and Bessel
  models. These covariance functions are defined as follows:

  \describe{
    \item{Whittle-Matern}{\eqn{\rho(h) = sill
	\frac{2^{1-smooth}}{\Gamma(smooth)} \left(\frac{h}{range}
	\right)^{smooth} K_{smooth}\left(\frac{h}{range} \right)}{rho(h) =
	sill * 2^(1-smooth) / \Gamma(smooth) * (h/range)^smooth * K_smooth
	(h / range)}}
    \item{Powered Exponential}{\eqn{\rho(h) = sill \exp \left[-
	\left(\frac{h}{range} \right)^{smooth} \right]}{rho(h) = sill *
	exp[-(h/range)^smooth]}}
    \item{Cauchy}{\eqn{\rho(h) = sill \left[1 + \left(\frac{h}{range}
	\right)^2 \right]^{-smooth}}{rho(h) = sill * [1 +
	(h/range)^2]^(-smooth)}}
    \item{Bessel}{\eqn{\rho(h) = sill \left(\frac{2
    range}{h}\right)^{smooth} \Gamma(smooth + 1)
    J_{smooth}\left(\frac{h}{range} \right)}{rho(h) = sill (2 range /
    h)^(smooth) \Gamma(smooth + 1) J_smooth(h/range)}}
  }
  where \eqn{\Gamma}{\Gamma} is the gamma function,
    \eqn{K_{smooth}}{K_smooth} is the modified Bessel function of order
    \eqn{smooth}{smooth} and \eqn{J_{smooth}}{J_smooth} is the Bessel
    function of order \eqn{smooth}{smooth}.
}
\examples{
## 1- Calling covariance using fixed covariance parameters
covariance(sill = 1, range = 1, smooth = 0.5, cov.mod = "whitmat")
covariance(sill = 0.5, range = 1, smooth = 0.5, cov.mod = "whitmat",
  dist = seq(0,5, 0.2), plot = FALSE)

## 2- Calling covariance from a fitted model
##Define the coordinate of each location
n.site <- 30
locations <- matrix(runif(2*n.site, 0, 10), ncol = 2)
colnames(locations) <- c("lon", "lat")

##Simulate a max-stable process - with unit Frechet margins
data <- rmaxstab(30, locations, cov.mod = "whitmat", sill = 1, range =
3, smooth = 0.5) 

##Fit a max-stable model
fitted <- fitmaxstab(data, locations, "whitmat", std.err.type = "none")
covariance(fitted, ylim = c(0, 1))
covariance(sill = 1, range = 3, smooth = .5, cov.mod = "whitmat", add =
TRUE, col = 3)
title("Whittle-Matern covariance function")
legend("topright", c("Theo.", "Fitted"), lty = 1, col = c(3,1), inset =
.05)
}
\author{Mathieu Ribatet}
\keyword{htest}