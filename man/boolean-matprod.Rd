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
\alias{\%&\%,ANY,sparseVector-method}
\alias{\%&\%,ANY,vector-method}
\alias{\%&\%,CsparseMatrix,CsparseMatrix-method}
\alias{\%&\%,CsparseMatrix,RsparseMatrix-method}
\alias{\%&\%,CsparseMatrix,TsparseMatrix-method}
\alias{\%&\%,CsparseMatrix,denseMatrix-method}
\alias{\%&\%,CsparseMatrix,diagonalMatrix-method}
\alias{\%&\%,CsparseMatrix,matrix-method}
\alias{\%&\%,CsparseMatrix,vector-method}
\alias{\%&\%,Matrix,ANY-method}
\alias{\%&\%,Matrix,indMatrix-method}
\alias{\%&\%,Matrix,pMatrix-method}
\alias{\%&\%,Matrix,sparseVector-method}
\alias{\%&\%,RsparseMatrix,CsparseMatrix-method}
\alias{\%&\%,RsparseMatrix,RsparseMatrix-method}
\alias{\%&\%,RsparseMatrix,TsparseMatrix-method}
\alias{\%&\%,RsparseMatrix,denseMatrix-method}
\alias{\%&\%,RsparseMatrix,diagonalMatrix-method}
\alias{\%&\%,RsparseMatrix,matrix-method}
\alias{\%&\%,RsparseMatrix,vector-method}
\alias{\%&\%,TsparseMatrix,CsparseMatrix-method}
\alias{\%&\%,TsparseMatrix,RsparseMatrix-method}
\alias{\%&\%,TsparseMatrix,TsparseMatrix-method}
\alias{\%&\%,TsparseMatrix,denseMatrix-method}
\alias{\%&\%,TsparseMatrix,diagonalMatrix-method}
\alias{\%&\%,TsparseMatrix,matrix-method}
\alias{\%&\%,TsparseMatrix,vector-method}
\alias{\%&\%,denseMatrix,CsparseMatrix-method}
\alias{\%&\%,denseMatrix,RsparseMatrix-method}
\alias{\%&\%,denseMatrix,TsparseMatrix-method}
\alias{\%&\%,denseMatrix,denseMatrix-method}
\alias{\%&\%,denseMatrix,diagonalMatrix-method}
\alias{\%&\%,denseMatrix,matrix-method}
\alias{\%&\%,denseMatrix,vector-method}
\alias{\%&\%,diagonalMatrix,CsparseMatrix-method}
\alias{\%&\%,diagonalMatrix,RsparseMatrix-method}
\alias{\%&\%,diagonalMatrix,TsparseMatrix-method}
\alias{\%&\%,diagonalMatrix,denseMatrix-method}
\alias{\%&\%,diagonalMatrix,diagonalMatrix-method}
\alias{\%&\%,diagonalMatrix,matrix-method}
\alias{\%&\%,diagonalMatrix,vector-method}
\alias{\%&\%,indMatrix,Matrix-method}
\alias{\%&\%,indMatrix,indMatrix-method}
\alias{\%&\%,indMatrix,matrix-method}
\alias{\%&\%,indMatrix,pMatrix-method}
\alias{\%&\%,indMatrix,vector-method}
\alias{\%&\%,matrix,ANY-method}
\alias{\%&\%,matrix,CsparseMatrix-method}
\alias{\%&\%,matrix,RsparseMatrix-method}
\alias{\%&\%,matrix,TsparseMatrix-method}
\alias{\%&\%,matrix,denseMatrix-method}
\alias{\%&\%,matrix,diagonalMatrix-method}
\alias{\%&\%,matrix,indMatrix-method}
\alias{\%&\%,matrix,matrix-method}
\alias{\%&\%,matrix,pMatrix-method}
\alias{\%&\%,matrix,sparseVector-method}
\alias{\%&\%,matrix,vector-method}
\alias{\%&\%,pMatrix,Matrix-method}
\alias{\%&\%,pMatrix,indMatrix-method}
\alias{\%&\%,pMatrix,matrix-method}
\alias{\%&\%,pMatrix,pMatrix-method}
\alias{\%&\%,pMatrix,vector-method}
\alias{\%&\%,sparseVector,ANY-method}
\alias{\%&\%,sparseVector,Matrix-method}
\alias{\%&\%,sparseVector,matrix-method}
\alias{\%&\%,sparseVector,sparseVector-method}
\alias{\%&\%,sparseVector,vector-method}
\alias{\%&\%,vector,ANY-method}
\alias{\%&\%,vector,CsparseMatrix-method}
\alias{\%&\%,vector,RsparseMatrix-method}
\alias{\%&\%,vector,TsparseMatrix-method}
\alias{\%&\%,vector,denseMatrix-method}
\alias{\%&\%,vector,diagonalMatrix-method}
\alias{\%&\%,vector,indMatrix-method}
\alias{\%&\%,vector,matrix-method}
\alias{\%&\%,vector,pMatrix-method}
\alias{\%&\%,vector,sparseVector-method}
\alias{\%&\%,vector,vector-method}
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
    \item{\code{signature(x = "nMatrix", y = "nMatrix")}}{ }
    \item{\code{signature(x = "nMatrix", y = "nsparseMatrix")}}{ }
    \item{\code{signature(x = "nsparseMatrix", y = "nMatrix")}}{ }
    \item{\code{signature(x = "nsparseMatrix", y = "nsparseMatrix")}}{ }
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
