\name{pack-methods}
\title{Representation of Packed and Unpacked Dense Matrices}
%
\docType{methods}
\keyword{array}
\keyword{methods}
%
\alias{pack}
\alias{pack-methods}
\alias{unpack}
\alias{unpack-methods}
%
\alias{pack,dgeMatrix-method}
\alias{pack,lgeMatrix-method}
\alias{pack,matrix-method}
\alias{pack,ngeMatrix-method}
\alias{pack,packedMatrix-method}
\alias{pack,sparseMatrix-method}
\alias{pack,unpackedMatrix-method}
%
\alias{unpack,matrix-method}
\alias{unpack,packedMatrix-method}
\alias{unpack,sparseMatrix-method}
\alias{unpack,unpackedMatrix-method}
%
\description{
  \code{pack()} coerces dense symmetric and dense triangular matrices
  from unpacked format (storing the full matrix) to packed format
  (storing only one of the upper and lower triangles).  \code{unpack()}
  performs the reverse coercion. The two formats are formalized
  by the virtual classes \code{"\linkS4class{packedMatrix}"} and
  \code{"\linkS4class{unpackedMatrix}"}.
}
\usage{
pack(x, \dots)
\S4method{pack}{dgeMatrix}(x, symmetric = NA, upperTri = NA, \dots)
\S4method{pack}{lgeMatrix}(x, symmetric = NA, upperTri = NA, \dots)
\S4method{pack}{ngeMatrix}(x, symmetric = NA, upperTri = NA, \dots)
\S4method{pack}{matrix}(x, symmetric = NA, upperTri = NA, \dots)

unpack(x, \dots)
}
\arguments{
  \item{x}{A dense symmetric or dense triangular matrix.
    \describe{
      \item{For \code{pack()}:}{typically an \code{"unpackedMatrix"}
	or a standard \code{"matrix"}, though \code{"packedMatrix"}
	are allowed and returned unchanged.}
      \item{For \code{unpack()}:}{typically a \code{"packedMatrix"},
	though \code{"unpackedMatrix"} are allowed and returned unchanged.}
    }
  }
  \item{symmetric}{logical (including \code{NA}) optionally
    indicating whether \code{x} is symmetric (or triangular).}
  \item{upperTri}{(for triangular \code{x} only) logical
    (including \code{NA}) indicating whether \code{x} is
    upper (or lower) triangular.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\value{
  \describe{
    \item{For \code{pack()}:}{a \code{"packedMatrix"} giving
      the condensed representation of \code{x}.}
    \item{For \code{unpack()}:}{an \code{"unpackedMatrix"} giving
      the full storage representation of \code{x}.}
  }
}
\details{
  \code{pack(x)} checks matrices \code{x} \emph{not} inheriting from
  one of the virtual classes \code{"\linkS4class{symmetricMatrix}"}
  \code{"\linkS4class{triangularMatrix}"} for symmetry
  (via \code{\link[=isSymmetric-methods]{isSymmetric}()})
  then for upper and lower triangularity
  (via \code{\link{isTriangular}()}) in order to identify a suitable
  coercion.  Setting one or both of \code{symmetric} and \code{upperTri}
  to \code{TRUE} or \code{FALSE} rather than \code{NA} allows skipping
  of irrelevant tests for large matrices known to be symmetric or
  (upper or lower) triangular.

  Users should \emph{not} assume that \code{pack()} and \code{unpack()}
  are inverse operations.  Specifically, \code{y <- unpack(pack(x))}
  may not reproduce an \code{"unpackedMatrix"} \code{x} in the sense of
  \code{\link{identical}()}.  See the examples.
}
\examples{
\dontshow{ % for R_DEFAULT_PACKAGES=NULL
library(utils, pos = "package:base", verbose = FALSE)
}
showMethods("pack")
(s <- crossprod(matrix(sample(15), 5,3))) # traditional symmetric matrix
(sp <- pack(s))
mt <- as.matrix(tt <- tril(s))
(pt <- pack(mt))
stopifnot(identical(pt, pack(tt)),
	  dim(s ) == dim(sp), all(s  == sp),
	  dim(mt) == dim(pt), all(mt == pt), all(mt == tt))

showMethods("unpack")
(cp4 <- chol(Hilbert(4))) # is triangular
tp4 <- pack(cp4) # [t]riangular [p]acked
str(tp4)
(unpack(tp4))
stopifnot(identical(tp4, pack(unpack(tp4))))

z1 <- new("dsyMatrix", Dim = c(2L, 2L), x = as.double(1:4), uplo = "U")
z2 <- unpack(pack(z1))
stopifnot(!identical(z1, z2), # _not_ identical
          all(z1 == z2)) # but mathematically equal
cbind(z1@x, z2@x) # (unused!) lower triangle is "lost" in translation
}
