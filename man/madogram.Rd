\name{madogram}
\alias{madogram}
\title{Computes madograms}
\description{Computes the madogram for max-stable processes.}
\usage{
madogram(data, coord, bins, gev.param = c(0, 1, 0), which =
c("mado", "ext"), xlab, ylab, angles = NULL, \dots)
}
\arguments{
  \item{data}{A matrix representing the data. Each column corresponds to
    one location.}
  \item{coord}{A matrix that gives the coordinates of each
    location. Each row corresponds to one location.}
  \item{bins}{Numeric vector giving the bins for which the madogram is
    computed.}
  \item{gev.param}{Numeric vector of length 3 specifying the location,
    scale and shape parameters for the GEV.}
  \item{which}{A character vector of maximum size 2. It specifies if the
    madogram and/or the extremal coefficient functions have to be
    plotted.}
  \item{xlab,ylab}{The x-axis and y-axis labels. May be missing. Note
    that \code{ylab} must have the same length has \code{which}.}
  \item{angles}{A numeric vector. A partition of the interval
    \eqn{(-\pi, \pi)}{(-\pi, \pi)} to help detecting anisotropy.}
  \item{\dots}{Additional options to be passed to the \code{\link{plot}}
    function.} 
}
\value{
  A graphic and (invisibly) a matrix with the lag distances, the
  estimated madogram values and extremal coefficients. 
}
\details{
  Let \eqn{Z(x)}{Z(x)} be a spatial and stationary process. The madogram
  is defined as follows:

  \deqn{\nu(h) = \frac{1}{2}\mbox{E}\left[|Z(x+h) - Z(x)|
    \right]}{\nu(h) = 0.5 * E[|Z(x+h) - Z(x)|]}

  If now \eqn{Z(x)}{Z(x)} is a stationary max-stable spatial random
  field with GEV marginals. Provided the GEV shape parameter
  \eqn{\xi}{\xi} is such that \eqn{\xi < 1}{\xi <1}. The extremal
  coefficient \eqn{\theta(h)}{\theta(h)} satisfies:

  \deqn{\theta(h) = \begin{cases} u_\beta \left(\mu +
    \frac{\nu(h)}{\Gamma(1 - \xi)} \right), & \xi \neq 0\\
    \exp\left(\frac{\nu(h)}{\sigma}\right), & \xi = 0 \end{cases}}{u_\beta (\mu +
    \nu(h) / \Gamma(1 - \xi)), if \xi < 1, exp(\nu(h)/\sigma), otherwise}
  where \eqn{\Gamma(\cdot)}{\Gamma(.)} is the gamma function and
  \eqn{u_\beta}{u_\beta} is defined as follows:

  \deqn{u_\beta(u) = \left(1 + \xi \frac{u - \mu}{\sigma}
    \right)_+^{1/\xi}}{(1 + \xi (u - \mu) / \sigma )_+^{1/\xi}} 
  and \eqn{\beta = (\mu, \sigma, \xi)}{\beta= (\mu, \sigma, \xi)} i.e
    the vector of the GEV parameters.
}
\references{
  Cooley, D., Naveau, P. and Poncet, P. (2006) Variograms for
  spatial max-stable random fields. \emph{Dependence in Probability and
    Statistics}, 373--390.
}
\seealso{
  \code{\link{fitextcoeff}}
}
\examples{
require(RandomFields)
n.site <- 50
locations <- matrix(runif(2*n.site, 0, 10), ncol = 2)
colnames(locations) <- c("lon", "lat")

##Simulate a max-stable process - with unit Frechet margins
ms0 <- MaxStableRF(locations[,1], locations[,2], grid=FALSE, model="wh",
                   param=c(0,1,0,1, 2), maxstable="extr",
                   n = 40)
ms0 <- t(ms0)

##Compute the madogram
madogram(ms0, locations, angles = seq(-pi,pi, pi/4))
}
\author{Mathieu Ribatet}
\keyword{hplot}
