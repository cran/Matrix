\name{USCounties}
\title{Contiguity Matrix of U.S. Counties}
%
\docType{data}
\keyword{datasets}
%
\alias{USCounties}
%
\description{
  This matrix gives the contiguities of 3111 U.S. counties,
  using the queen criterion of at least one shared vertex
  or edge.
}
\usage{data(USCounties)}
\format{
  A \eqn{3111 \times 3111}{3111-by-3111} sparse, symmetric
  matrix of class \code{\linkS4class{dsCMatrix}}, with 9101
  nonzero entries.
}
\source{
  GAL lattice file \file{usc_q.GAL}
  (retrieved in 2008 from
  \file{http://sal.uiuc.edu/weights/zips/usc.zip}
  with permission from Luc Anselin for use and distribution)
  was read into \R{} using function \code{read.gal}
  from package \CRANpkg{spdep}.

  Neighbour lists were augmented with row-standardized
  (and then symmetrized) spatial weights, using functions
  \code{nb2listw} and \code{similar.listw} from packages
  \CRANpkg{spdep} and \CRANpkg{spatialreg}.
  The resulting \code{listw} object was coerced to class
  \code{\linkS4class{dsTMatrix}}
  using \code{as_dsTMatrix_listw} from \CRANpkg{spatialreg},
  and subsequently to class \code{\linkS4class{dsCMatrix}}.
}
\references{
  Ord, J. K. (1975).
  Estimation methods for models of spatial interaction.
  \emph{Journal of the American Statistical Association},
  \emph{70}(349), 120-126.
  \doi{10.2307/2285387}
}
\examples{
\dontshow{ % for R_DEFAULT_PACKAGES=NULL
library(stats, pos = "package:base", verbose = FALSE)
library(utils, pos = "package:base", verbose = FALSE)
}
data(USCounties, package = "Matrix")
(n <- ncol(USCounties))
I <- .symDiagonal(n)

set.seed(1)
r <- 50L
rho <- 1 / runif(r, 0, 0.5)

system.time(MJ0 <- sapply(rho, function(mult)
    determinant(USCounties + mult * I, logarithm = TRUE)$modulus))

## Can be done faster by updating the Cholesky factor:

C1 <- Cholesky(USCounties, Imult = 2)
system.time(MJ1 <- sapply(rho, function(mult)
    determinant(update(C1, USCounties, mult), sqrt = FALSE)$modulus))
stopifnot(all.equal(MJ0, MJ1))

C2 <- Cholesky(USCounties, super = TRUE, Imult = 2)
system.time(MJ2 <- sapply(rho, function(mult)
    determinant(update(C2, USCounties, mult), sqrt = FALSE)$modulus))
stopifnot(all.equal(MJ0, MJ2))
}
