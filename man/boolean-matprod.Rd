\name{boolmatmult-methods}
\title{Boolean Arithmetic Matrix Products: \code{\%&\%} and Methods}
%
\docType{methods}
\keyword{algebra}
\keyword{array}
\keyword{logic}
\keyword{methods}
%
\alias{\%&\%}
\alias{\%&\%-methods}
\alias{boolmatmult-methods}
%
\alias{\%&\%,ANY,ANY-method}
\alias{\%&\%,ANY,Matrix-method}
\alias{\%&\%,ANY,matrix-method}
\alias{\%&\%,CsparseMatrix,RsparseMatrix-method}
\alias{\%&\%,CsparseMatrix,TsparseMatrix-method}
\alias{\%&\%,CsparseMatrix,diagonalMatrix-method}
\alias{\%&\%,CsparseMatrix,mMatrix-method}
\alias{\%&\%,Matrix,ANY-method}
\alias{\%&\%,Matrix,Matrix-method}
\alias{\%&\%,Matrix,indMatrix-method}
\alias{\%&\%,Matrix,pMatrix-method}
\alias{\%&\%,RsparseMatrix,CsparseMatrix-method}
\alias{\%&\%,RsparseMatrix,RsparseMatrix-method}
\alias{\%&\%,RsparseMatrix,TsparseMatrix-method}
\alias{\%&\%,RsparseMatrix,diagonalMatrix-method}
\alias{\%&\%,RsparseMatrix,mMatrix-method}
\alias{\%&\%,TsparseMatrix,CsparseMatrix-method}
\alias{\%&\%,TsparseMatrix,RsparseMatrix-method}
\alias{\%&\%,TsparseMatrix,TsparseMatrix-method}
\alias{\%&\%,TsparseMatrix,diagonalMatrix-method}
\alias{\%&\%,TsparseMatrix,mMatrix-method}
\alias{\%&\%,denseMatrix,denseMatrix-method}
\alias{\%&\%,denseMatrix,diagonalMatrix-method}
\alias{\%&\%,diagonalMatrix,CsparseMatrix-method}
\alias{\%&\%,diagonalMatrix,RsparseMatrix-method}
\alias{\%&\%,diagonalMatrix,TsparseMatrix-method}
\alias{\%&\%,diagonalMatrix,denseMatrix-method}
\alias{\%&\%,diagonalMatrix,diagonalMatrix-method}
\alias{\%&\%,diagonalMatrix,matrix-method}
\alias{\%&\%,indMatrix,Matrix-method}
\alias{\%&\%,indMatrix,indMatrix-method}
\alias{\%&\%,indMatrix,matrix-method}
\alias{\%&\%,indMatrix,pMatrix-method}
\alias{\%&\%,mMatrix,CsparseMatrix-method}
\alias{\%&\%,mMatrix,RsparseMatrix-method}
\alias{\%&\%,mMatrix,TsparseMatrix-method}
\alias{\%&\%,mMatrix,sparseMatrix-method}
\alias{\%&\%,mMatrix,sparseVector-method}
\alias{\%&\%,matrix,ANY-method}
\alias{\%&\%,matrix,diagonalMatrix-method}
\alias{\%&\%,matrix,indMatrix-method}
\alias{\%&\%,matrix,matrix-method}
\alias{\%&\%,matrix,pMatrix-method}
\alias{\%&\%,nCsparseMatrix,nCsparseMatrix-method}
\alias{\%&\%,nCsparseMatrix,nsparseMatrix-method}
\alias{\%&\%,nMatrix,nMatrix-method}
\alias{\%&\%,nMatrix,nsparseMatrix-method}
\alias{\%&\%,nsparseMatrix,nCsparseMatrix-method}
\alias{\%&\%,nsparseMatrix,nMatrix-method}
\alias{\%&\%,nsparseMatrix,nsparseMatrix-method}
\alias{\%&\%,numLike,sparseVector-method}
\alias{\%&\%,pMatrix,Matrix-method}
\alias{\%&\%,pMatrix,indMatrix-method}
\alias{\%&\%,pMatrix,matrix-method}
\alias{\%&\%,pMatrix,pMatrix-method}
\alias{\%&\%,sparseMatrix,mMatrix-method}
\alias{\%&\%,sparseMatrix,sparseMatrix-method}
\alias{\%&\%,sparseVector,mMatrix-method}
\alias{\%&\%,sparseVector,numLike-method}
\alias{\%&\%,sparseVector,sparseVector-method}
%
\description{
  For boolean or \dQuote{patter\bold{n}} matrices, i.e., \R objects of
  class \code{\linkS4class{nMatrix}}, it is natural to allow matrix
  products using boolean instead of numerical arithmetic.

  In package \pkg{Matrix}, we use the binary operator \code{\%&\%} (aka
  \dQuote{infix}) function) for this and provide methods for all our
  matrices and the traditional \R matrices (see \code{\link{matrix}}).
}
\section{Methods}{
  We provide methods for both the \dQuote{traditional} (\R base) matrices
  and numeric vectors and conceptually all matrices and
  \code{\linkS4class{sparseVector}}s in package \pkg{Matrix}.
  \describe{
    \item{\code{signature(x = "ANY", y = "ANY")}}{ }
    \item{\code{signature(x = "ANY", y = "Matrix")}}{ }
    \item{\code{signature(x = "Matrix", y = "ANY")}}{ }
    \item{\code{signature(x = "mMatrix", y = "mMatrix")}}{ }
    \item{\code{signature(x = "nMatrix", y = "nMatrix")}}{ }
    \item{\code{signature(x = "nMatrix", y = "nsparseMatrix")}}{ }
    \item{\code{signature(x = "nsparseMatrix", y = "nMatrix")}}{ }
    \item{\code{signature(x = "nsparseMatrix", y = "nsparseMatrix")}}{ }
    \item{\code{signature(x = "sparseVector", y = "mMatrix")}}{ }
    \item{\code{signature(x = "mMatrix", y = "sparseVector")}}{ }
    \item{\code{signature(x = "sparseVector", y = "sparseVector")}}{ }
  }% {describe}
}% {Methods}
\note{
  These boolean arithmetic matrix products had been newly
  introduced for \pkg{Matrix} 1.2.0 (March 2015).  Its implementation
  has still not been tested extensively.

  Originally, it was left unspecified how non-structural zeros, i.e., \code{0}'s
  as part of the \code{M@x} slot should be treated for numeric
  (\code{"\linkS4class{dMatrix}"}) and logical (\code{"\linkS4class{lMatrix}"})
  sparse matrices.  We now specify that boolean matrix products should behave as if
  applied to \code{\link{drop0}(M)}, i.e., as if dropping such zeros from
  the matrix before using it.
  \cr
  Equivalently, for all matrices \code{M}, boolean arithmetic should work as if
  applied to \code{M != 0} (or \code{M != FALSE}).

  The current implementation ends up coercing both \code{x} and \code{y} to
  (virtual) class \code{\linkS4class{nsparseMatrix}} which may be quite inefficient
  for dense matrices.  A future implementation may well return a matrix
  with \bold{different} class, but the \dQuote{same} content, i.e., the
  same matrix entries \eqn{m_ij}{m[i,j]}.
}
\seealso{
  \code{\link{\%*\%}}, \code{crossprod()}, or \code{\link{tcrossprod}()},
  for (regular) matrix product methods.
}
\value{
  a pattern matrix, i.e., inheriting from \code{"\linkS4class{nMatrix}"},
  or an \code{"\linkS4class{ldiMatrix}"} in case of a diagonal matrix.
}
\examples{
\dontshow{ % for R_DEFAULT_PACKAGES=NULL
library(stats, pos = "package:base", verbose = FALSE)
}
set.seed(7)
L <- Matrix(rnorm(20) > 1,    4,5)
(N <- as(L, "nMatrix"))
L. <- L; L.[1:2,1] <- TRUE; L.@x[1:2] <- FALSE; L. # has "zeros" to drop0()
D <- Matrix(round(rnorm(30)), 5,6) # -> values in -1:1 (for this seed)
L \%&\% D
stopifnot(identical(L \%&\% D, N \%&\% D),
          all(L \%&\% D == as((L \%*\% abs(D)) > 0, "sparseMatrix")))

## cross products , possibly with  boolArith = TRUE :
crossprod(N)     # -> sparse patter'n' (TRUE/FALSE : boolean arithmetic)
crossprod(N  +0) # -> numeric Matrix (with same "pattern")
stopifnot(all(crossprod(N) == t(N) \%&\% N),
          identical(crossprod(N), crossprod(N +0, boolArith=TRUE)),
          identical(crossprod(L), crossprod(N   , boolArith=FALSE)))
crossprod(D, boolArith =  TRUE) # pattern: "nsCMatrix"
crossprod(L, boolArith =  TRUE) #  ditto
crossprod(L, boolArith = FALSE) # numeric: "dsCMatrix"
}
