\name{kronecker-methods}
\title{Methods for Function 'kronecker()' in Package 'Matrix'}
%
\docType{methods}
\keyword{algebra}
\keyword{arith}
\keyword{array}
\keyword{methods}
%
\alias{kronecker}
\alias{kronecker-methods}
%
\alias{kronecker,CsparseMatrix,CsparseMatrix-method}
\alias{kronecker,CsparseMatrix,Matrix-method}
\alias{kronecker,CsparseMatrix,diagonalMatrix-method}
\alias{kronecker,Matrix,matrix-method}
\alias{kronecker,Matrix,vector-method}
\alias{kronecker,RsparseMatrix,Matrix-method}
\alias{kronecker,RsparseMatrix,RsparseMatrix-method}
\alias{kronecker,RsparseMatrix,diagonalMatrix-method}
\alias{kronecker,TsparseMatrix,Matrix-method}
\alias{kronecker,TsparseMatrix,TsparseMatrix-method}
\alias{kronecker,TsparseMatrix,diagonalMatrix-method}
\alias{kronecker,denseMatrix,Matrix-method}
\alias{kronecker,denseMatrix,denseMatrix-method}
\alias{kronecker,diagonalMatrix,CsparseMatrix-method}
\alias{kronecker,diagonalMatrix,Matrix-method}
\alias{kronecker,diagonalMatrix,RsparseMatrix-method}
\alias{kronecker,diagonalMatrix,TsparseMatrix-method}
\alias{kronecker,diagonalMatrix,diagonalMatrix-method}
\alias{kronecker,diagonalMatrix,indMatrix-method}
\alias{kronecker,indMatrix,Matrix-method}
\alias{kronecker,indMatrix,diagonalMatrix-method}
\alias{kronecker,indMatrix,indMatrix-method}
\alias{kronecker,matrix,Matrix-method}
\alias{kronecker,vector,Matrix-method}
%
\description{
  Computes Kronecker products for objects inheriting from
  \code{"\linkS4class{Matrix}"}.

  In order to preserver sparseness, we treat \code{0 * NA} as \code{0},
  not as \code{\link{NA}} as usually in \R (and as used for the
  \pkg{base} function \code{\link[base]{kronecker}}).
}
\section{Methods}{
  \describe{

    \item{kronecker}{\code{signature(X = "Matrix", Y = "ANY")} .......}
    \item{kronecker}{\code{signature(X = "ANY", Y = "Matrix")} .......}
    \item{kronecker}{\code{signature(X = "diagonalMatrix", Y = "ANY")} .......}
    \item{kronecker}{\code{signature(X = "sparseMatrix", Y = "ANY")} .......}
    \item{kronecker}{\code{signature(X = "TsparseMatrix", Y = "TsparseMatrix")} .......}
    \item{kronecker}{\code{signature(X = "dgTMatrix", Y = "dgTMatrix")} .......}
    \item{kronecker}{\code{signature(X = "dtTMatrix", Y = "dtTMatrix")} .......}
    \item{kronecker}{\code{signature(X = "indMatrix", Y = "indMatrix")} .......}

  }
}
\examples{
(t1 <- spMatrix(5,4, x= c(3,2,-7,11), i= 1:4, j=4:1)) #  5 x  4
(t2 <- kronecker(Diagonal(3, 2:4), t1))               # 15 x 12

## should also work with special-cased logical matrices
l3 <- upper.tri(matrix(,3,3))
M <- Matrix(l3)
(N <- as(M, "nsparseMatrix")) # "ntCMatrix" (upper triangular)
N2 <- as(N, "generalMatrix")  # (lost "t"riangularity)
MM <- kronecker(M,M)
NN <- kronecker(N,N) # "dtTMatrix" i.e. did keep
NN2 <- kronecker(N2,N2)
stopifnot(identical(NN,MM),
          is(NN2, "sparseMatrix"), all(NN2 == NN),
          is(NN, "triangularMatrix"))
}
\keyword{methods}
\keyword{array}

