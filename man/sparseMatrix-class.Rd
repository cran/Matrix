\name{sparseMatrix-class}
\docType{class}
\title{Virtual Class "sparseMatrix" --- Mother of Sparse Matrices}
\alias{sparseMatrix-class}
%% Group methods
\alias{Compare,sparseMatrix,sparseMatrix-method}
\alias{Arith,sparseMatrix,sparseMatrix-method}
\alias{Arith,sparseMatrix,numeric-method}
\alias{Arith,numeric,sparseMatrix-method}
\alias{Math,sparseMatrix-method}
%\alias{Math2,sparseMatrix,numeric-method}
\alias{Ops,diagonalMatrix,sparseMatrix-method}
\alias{Ops,sparseMatrix,diagonalMatrix-method}
%
\alias{cbind2,sparseMatrix,sparseMatrix-method}
\alias{cbind2,sparseMatrix,matrix-method}
\alias{cbind2,matrix,sparseMatrix-method}
\alias{cbind2,sparseMatrix,numeric-method}
\alias{cbind2,numeric,sparseMatrix-method}
\alias{rbind2,sparseMatrix,sparseMatrix-method}
\alias{rbind2,sparseMatrix,matrix-method}
\alias{rbind2,matrix,sparseMatrix-method}
\alias{rbind2,sparseMatrix,numeric-method}
\alias{rbind2,numeric,sparseMatrix-method}
%
\alias{coerce,ANY,sparseMatrix-method}
\alias{coerce,sparseMatrix,generalMatrix-method}
\alias{-,sparseMatrix,missing-method}
\alias{diag,sparseMatrix-method}
\alias{show,sparseMatrix-method}
\alias{colMeans,sparseMatrix-method}
\alias{colSums,sparseMatrix-method}
\alias{rowMeans,sparseMatrix-method}
\alias{rowSums,sparseMatrix-method}
% graph stuff
\alias{coerce,graph,sparseMatrix-method}
\alias{coerce,graph,CsparseMatrix-method}
\alias{coerce,graphAM,sparseMatrix-method}
\alias{coerce,graphNEL,CsparseMatrix-method}
\alias{coerce,graphNEL,TsparseMatrix-method}
\alias{coerce,sparseMatrix,graph-method}
\alias{coerce,sparseMatrix,graphNEL-method}
\alias{coerce,TsparseMatrix,graphNEL-method}
\alias{coerce,dgTMatrix,graphNEL-method}
%
\description{Virtual Mother Class of All Sparse Matrices}
\section{Slots}{
  \describe{
    \item{\code{Dim}:}{Object of class \code{"integer"} - the dimensions
     of the matrix - must be an integer vector with exactly two
     non-negative values.}
    \item{\code{Dimnames}:}{a list of length two - inherited from class
      \code{Matrix}, see \code{\linkS4class{Matrix}}.}
    \item{\code{factors}:}{Object of class \code{"list"} - a list
      of factorizations of the matrix. }
  }
}
\section{Extends}{
  Class \code{"Matrix"}, directly.
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "sparseMatrix")}: The
      \code{\link{show}} method for sparse matrices prints
      \emph{\dQuote{structural}} zeroes as \code{"."} using the
      non-exported \code{prSpMatrix} function which allows further
      customization such as replacing \code{"."} by \code{" "} (blank).

      Note that \code{\link{options}(max.print)} will influence how many
      entries of large matrices are printed at all.}
  }
}
\note{
  In method selection for multiplication operations (i.e. \code{\%*\%}
  and the two-argument form of \code{\link[base]{crossprod}})
  the sparseMatrix class takes precedence in the sense that if one
  operand is a sparse matrix and the other is any type of dense matrix
  then the dense matrix is coerced to a \code{dgeMatrix} and the
  appropriate sparse matrix method is used.
}
%\author{Martin}
\examples{
showClass("sparseMatrix") ## and look at the help() of its subclasses
M <- Matrix(0, 10000, 100)
M[1,1] <- M[2,3] <- 3.14
M
\dontshow{% show() was slow in 0.9975-8 because of slow adjust="fancy"
op <- options(max.print = 25000, width = 80)
sink(print(tempfile()))
stopifnot((st <- system.time(show(M)))[1] < 1.0) # only 0.09 on cmath-3
sink()
st
options(op)
}
}
\keyword{classes}
