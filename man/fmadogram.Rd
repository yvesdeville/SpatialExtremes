\name{fmadogram}
\alias{fmadogram}
\title{Computes the F-madogram}
\description{Computes the F-madogram for max-stable processes.}
\usage{
fmadogram(data, coord, n.bins, which = c("mado", "ext"), xlab, ylab,
angles = NULL, marge = "mle", ...)
}
\arguments{
  \item{data}{A matrix representing the data. Each column corresponds to
    one location.}
  \item{coord}{A matrix that gives the coordinates of each
    location. Each row corresponds to one location.}
  \item{n.bins}{The number of bins to be used. If missing, pairwise
    F-madogram estimates will be computed.}
  \item{which}{A character vector of maximum size 2. It specifies if the
    madogram and/or the extremal coefficient functions have to be
    plotted.}
  \item{xlab,ylab}{The x-axis and y-axis labels. May be missing. Note
    that \code{ylab} must have the same length has \code{which}.}
  \item{angles}{A numeric vector. A partition of the interval
    \eqn{(-\pi, \pi)}{(-\pi, \pi)} to help detecting anisotropy.}
  \item{marge}{Character string. If 'emp', the probabilities of non
    exceedances are estimated using the empirical CDF. If 'mle'
    (default), maximum likelihood estimates are used.}
  \item{\dots}{Additional options to be passed to the \code{\link{plot}}
    function.} 
}
\value{
  A graphic and (invisibly) a matrix with the lag distances, the
  F-madogram  and extremal coefficient estimates. 
}
\details{
  Let \eqn{Z(x)}{Z(x)} be a stationary process. The F-madogram is
  defined as follows:

  \deqn{\nu(h) = \frac{1}{2}\mbox{E}\left[|F(Z(x+h)) - F(Z(x))|
    \right]}{\nu(h) = 0.5 * E[|F(Z(x+h)) - F(Z(x))|]}

  The extremal coefficient \eqn{\theta(h)}{\theta(h)} satisfies:

  \deqn{\theta(h) = \frac{1 + 2 \nu(h)}{1 - 2 \nu(h)}}{(1 + 2 nu(h)) /
    (1 - 2 nu(h)}
}
\references{
  Cooley, D., Naveau, P. and Poncet, P. (2006) Variograms for
  spatial max-stable random fields. \emph{Dependence in Probability and
    Statistics}, 373--390.
}
\seealso{
  \code{\link{madogram}}, \code{\link{lmadogram}}
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

##Compute the F-madogram
fmadogram(ms0, locations)
}
\author{Mathieu Ribatet}
\keyword{hplot}
