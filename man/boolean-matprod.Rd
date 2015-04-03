\name{\%&\%-methods}
\title{Boolean Arithmetic Matrix Products: \code{\%&\%} and Methods}
\docType{methods}
\alias{\%&\%}
\alias{\%&\%-methods}
\alias{\%&\%,ANY,ANY-method}
\alias{\%&\%,ANY,matrix-method}
\alias{\%&\%,ANY,Matrix-method}
\alias{\%&\%,matrix,ANY-method}
\alias{\%&\%,Matrix,ANY-method}
\alias{\%&\%,matrix,matrix-method}
\alias{\%&\%,Matrix,Matrix-method}
\alias{\%&\%,mMatrix,mMatrix-method}
\alias{\%&\%,mMatrix,nMatrix-method}
\alias{\%&\%,nMatrix,mMatrix-method}
\alias{\%&\%,nMatrix,nMatrix-method}
\alias{\%&\%,nMatrix,nsparseMatrix-method}
\alias{\%&\%,nsparseMatrix,nMatrix-method}
\alias{\%&\%,nsparseMatrix,nsparseMatrix-method}
\alias{\%&\%,nCsparseMatrix,nCsparseMatrix-method}
\alias{\%&\%,nCsparseMatrix,nsparseMatrix-method}
\alias{\%&\%,nsparseMatrix,nCsparseMatrix-method}
\alias{\%&\%,diagonalMatrix,geMatrix-method}
\alias{\%&\%,geMatrix,diagonalMatrix-method}
\alias{\%&\%,CsparseMatrix,diagonalMatrix-method}
\alias{\%&\%,diagonalMatrix,CsparseMatrix-method}
\alias{\%&\%,diagonalMatrix,diagonalMatrix-method}
\alias{\%&\%,diagonalMatrix,sparseMatrix-method}
\alias{\%&\%,sparseMatrix,diagonalMatrix-method}
\alias{\%&\%,mMatrix,sparseVector-method}
\alias{\%&\%,numLike,sparseVector-method}
\alias{\%&\%,sparseVector,mMatrix-method}
\alias{\%&\%,sparseVector,numLike-method}
\alias{\%&\%,sparseVector,sparseVector-method}

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
  The current implementation ends up coercing both \code{x} and \code{y}
  to (virtual) class \code{\linkS4class{nsparseMatrix}} which may be
  quite inefficient.  A future implementation may well return a matrix
  with \bold{different} class, but the \dQuote{same} content, i.e., the
  same matrix entries \eqn{m_ij}{m[i,j]}.
}
\value{
  a pattern matrix, i.e., inheriting from \code{"\linkS4class{nMatrix}"},
  or an \code{"\linkS4class{ldiMatrix}"} in case of a diagonal matrix.
}
\examples{
set.seed(7)
L <- Matrix(rnorm(20) > 1,    4,5)
(N <- as(L, "nMatrix"))
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
\keyword{methods}
