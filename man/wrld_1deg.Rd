\name{wrld_1deg}
\title{World 1-degree grid contiguity matrix}
\alias{wrld_1deg}
\docType{data}
\description{
  This matrix represents the distance-based contiguities of 15260
  one-degree grid cells of land areas. The representation is as a
  row standardised spatial weights matrix transformed to a symmetric
  matrix (see Ord (1975), p. 125).
}
\usage{data(wrld_1deg)}
\format{
  A \eqn{15260 ^2} symmetric sparse matrix of class
  \code{\linkS4class{dsCMatrix}} with 55973 non-zero entries.
}
\details{
  The data were created into \R using the coordinates of a
  \sQuote{SpatialPixels} object containing approximately one-degree grid
  cells for land areas only (world excluding Antarctica), using package
  \CRANpkg{spdep}'s
  \code{\link[spdep]{dnearneigh}} with a cutoff distance of \code{sqrt(2)},
  and row-standardised and transformed to symmetry using
  \code{\link[spdep]{nb2listw}} and \code{\link[spdep]{similar.listw}}.
  This spatial weights object was converted to a \code{\linkS4class{dsTMatrix}}
  using \code{\link[spdep]{as_dsTMatrix_listw}} and then coerced
  (column-compressed).
}
\source{
  The shoreline data was read into \R using \code{\link[maptools]{Rgshhs}}
  from the GSHHS coarse shoreline database distributed with the
  \CRANpkg{maptools} package, omitting Antarctica.  A matching approximately
  one-degree grid was generated using
  \code{\link[maptools]{Sobj_SpatialGrid}}, and the grids on land were
  found using the appropriate \code{\link[sp]{over}}% was {overlay} till 2016-02
  method for the \sQuote{SpatialPolygons} and \sQuote{SpatialGrid} objects,
  yielding a \sQuote{SpatialPixels} one containing only the grid cells with
  centres on land.
}
\references{
  Ord, J. K. (1975)
  Estimation methods for models of spatial interaction;
  \emph{Journal of the American Statistical Association} \bold{70}, 120--126.
}
\examples{
data(wrld_1deg)
(n <- ncol(wrld_1deg))
IM <- .symDiagonal(n)
doExtras <- interactive() || nzchar(Sys.getenv("R_MATRIX_CHECK_EXTRA")) ||
    identical("true", unname(Sys.getenv("R_PKG_CHECKING_doExtras")))
nn <- if(doExtras) 20 else 3
set.seed(1)
rho <- runif(nn, 0, 1)
system.time(MJ <- sapply(rho,
                   function(x) determinant(IM - x * wrld_1deg,
                                           logarithm = TRUE)$modulus))
nWC <- -wrld_1deg
C1 <- Cholesky(nWC, Imult = 2)
## Note that det(<CHMfactor>) = det(L) = sqrt(det(A))
## ====> log det(A) = log( det(L)^2 ) = 2 * log det(L) :
system.time(MJ1 <- n * log(rho) +
   sapply(rho, function(x) c(2* determinant(update(C1, nWC, 1/x))$modulus))
)
stopifnot(all.equal(MJ, MJ1))
C2 <- Cholesky(nWC, super = TRUE, Imult = 2)
system.time(MJ2 <- n * log(rho) +
   sapply(rho, function(x) c(2* determinant(update(C2, nWC, 1/x))$modulus))
)
system.time(MJ3 <- n * log(rho) + Matrix:::ldetL2up(C1, nWC, 1/rho))
system.time(MJ4 <- n * log(rho) + Matrix:::ldetL2up(C2, nWC, 1/rho))
stopifnot(all.equal(MJ, MJ2),
          all.equal(MJ, MJ3),
          all.equal(MJ, MJ4))
}
\keyword{datasets}
