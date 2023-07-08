\name{wrld_1deg}
\title{Contiguity Matrix of World One-Degree Grid Cells}
%
\docType{data}
\keyword{datasets}
%
\alias{wrld_1deg}
%
\description{
  This matrix gives the contiguities of 15260 one-degree
  grid cells of world land areas, using a criterion based
  on the great-circle distance between centers.
}
\usage{data(wrld_1deg)}
\format{
  A \eqn{15260 \times 15260}{15260-by-15260} sparse, symmetric
  matrix of class \code{\linkS4class{dsCMatrix}}, with 55973
  nonzero entries.
}
\source{
  Shoreline data were read into \R{} from the GSHHS database
  using function \code{Rgshhs} from package \CRANpkg{maptools}.
  Antarctica was excluded.  An approximately one-degree grid
  was generated using function \code{Sobj_SpatialGrid}, also
  from \CRANpkg{maptools}.  Grid cells with centers on land
  were identified using the \code{over} method for classes
  \code{SpatialPolygons} and \code{SpatialGrid}, defined in
  package \CRANpkg{sp}.  Neighbours of these were identified
  by passing the resulting \code{SpatialPixels} object to
  function \code{dnearneigh} from package \CRANpkg{spdep},
  using as a cut-off a great-circle distance of \code{sqrt(2)}
  kilometers between centers.

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
data(wrld_1deg, package = "Matrix")
(n <- ncol(wrld_1deg))
I <- .symDiagonal(n)

doExtras <- interactive() || nzchar(Sys.getenv("R_MATRIX_CHECK_EXTRA"))
set.seed(1)
r <- if(doExtras) 20L else 3L
rho <- 1 / runif(r, 0, 0.5)

system.time(MJ0 <- sapply(rho, function(mult)
    determinant(wrld_1deg + mult * I, logarithm = TRUE)$modulus))

## Can be done faster by updating the Cholesky factor:

C1 <- Cholesky(wrld_1deg, Imult = 2)
system.time(MJ1 <- sapply(rho, function(mult)
    determinant(update(C1, wrld_1deg, mult), sqrt = FALSE)$modulus))
stopifnot(all.equal(MJ0, MJ1))

C2 <- Cholesky(wrld_1deg, super = TRUE, Imult = 2)
system.time(MJ2 <- sapply(rho, function(mult)
    determinant(update(C2, wrld_1deg, mult), sqrt = FALSE)$modulus))
stopifnot(all.equal(MJ0, MJ2))
}
