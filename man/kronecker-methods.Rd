\name{kronecker-methods}
\docType{methods}
\title{Methods for Function 'kronecker()' in Package 'Matrix'}
\alias{kronecker-methods}
%
\alias{kronecker,ANY,Matrix-method}
\alias{kronecker,ANY,diagonalMatrix-method}
\alias{kronecker,ANY,sparseMatrix-method}
\alias{kronecker,Matrix,ANY-method}
\alias{kronecker,diagonalMatrix,ANY-method}
\alias{kronecker,sparseMatrix,ANY-method}
%
\alias{kronecker,TsparseMatrix,TsparseMatrix-method}
\alias{kronecker,dgTMatrix,dgTMatrix-method}
\alias{kronecker,dgTMatrix,dtTMatrix-method}
\alias{kronecker,dtTMatrix,dgTMatrix-method}
\alias{kronecker,dtTMatrix,dtTMatrix-method}
\alias{kronecker,dsparseMatrix,dsparseMatrix-method}
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

    \item{kronecker}{\code{signature(X = "diagonalMatrix", Y = "ANY")} .......}
    \item{kronecker}{\code{signature(X = "diagonalMatrix", Y = "ANY")} .......}
    \item{kronecker}{\code{signature(X = "diagonalMatrix", Y = "ANY")} .......}

}}
\keyword{methods}
\keyword{array}

