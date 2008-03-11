\name{anova}
\alias{anova}
\alias{anova.maxstab}
\title{Anova Tables}
\description{
  Computes analysis of deviance for objects of class ''maxstab''
}
\usage{
\method{anova}{maxstab}(object, object2, method = "RJ", square = "chol",
...)
}
\arguments{
  \item{object, object2}{Two objects of class 'maxstab'.}
  \item{method}{Character string. Must be one of "CB" or "RJ" for the
    Chandler and Bate or the Rotnitzky and Jewell approaches
    respectively. See function
    \code{\link[=profile.maxstab]{profile}}.}
  \item{square}{The choice for the matrix square root. This is only
    useful for the 'CB' method. Must be one of 'chol' (Cholesky) or
    'svd' (Singular Value Decomposition).}
  \item{\dots}{Other options to be passed to the \code{\link{anova}}
    function.}
}
\details{
  As ''maxstab'' objects are fitted using pairwise likelihood, the model
  is misspecified. As a consequence, the likelihood ratio statistic is
  no longer \eqn{\chi^2}{\chi^2} distributed. To compute the anova
  table, we use the methodology proposed by Rotnitzky and Jewell to
  adjust the distribution of the likelihood ratio statistic.
}
\value{
This function returns an object of class anova. These objects
represent analysis-of-deviance tables. 
}
\references{
  Chandler, R. E. and Bate, S. (2007) Inference for clustered data using
  the independence loglikelihood \emph{Biometrika}, \bold{94},
  167--183.
  
  Rotnitzky, A. and Jewell, N. (1990) Hypothesis testing of regression
  parameters in semiparametric generalized linear models for cluster
  correlated data. \emph{Biometrika} \bold{77}, 485--97.
}
\author{Mathieu Ribatet}
\seealso{\code{\link{fitmaxstab}},
  \code{\link[=profile.maxstab]{profile}}, \code{\link{TIC}}}
\examples{
require(RandomFields)

##Define the coordinates of each location
n.site <- 30
locations <- matrix(rnorm(2*n.site, sd = sqrt(.2)), ncol = 2)
colnames(locations) <- c("lon", "lat")

##Simulate a max-stable process - with unit Frechet margins
sigma <- matrix(c(100, 25, 25, 220),ncol = 2)
sigma.inv <- solve(sigma)
sqrtCinv <- t(chol(sigma.inv))
model <- list(list(model = "gauss", var = 1, aniso = sqrtCinv / 2))
ms0 <- MaxStableRF(locations[,1], locations[,2], grid=FALSE, model=model,
                   maxstable = "Bool", n = 50)
ms0 <- t(ms0)

##Fit three nested models
M0 <- fitmaxstab(ms0, locations, "gauss", fit.marge = FALSE)
M1 <- fitmaxstab(ms0, locations, "gauss", fit.marge = FALSE, cov11 =
100)
M2 <- fitmaxstab(ms0, locations, "gauss", fit.marge = FALSE, cov11 =
100, cov22 = 500, method = "BFGS")
anova(M0, M1) ## <- M1 seems better
anova(M0, M2) ## <- M0 is better
anova(M1, M2) ## <- M1 is clearly better

##To use the Chandler and Bate methodology:
anova(M0, M1, method = "CB")
}
\keyword{htest}
