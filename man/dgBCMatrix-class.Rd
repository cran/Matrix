\name{dgBCMatrix-class}
\docType{class}
\alias{dgBCMatrix-class}
\alias{coerce,dgBCMatrix,dgeMatrix-method}
\alias{coerce,dgBCMatrix,dgCMatrix-method}
\alias{coerce,dgBCMatrix,dgTMatrix-method}
\alias{coerce,dgBCMatrix,matrix-method}
\title{Class "dgBCMatrix" Real Sparse Blocked Column Compressed Matrix}
\description{
  The \code{dgCMatrix} class is a class of sparse blocked matrices in
  the compressed, sparse, column-oriented format.  In this
  implementation the non-zero elements in the columns are sorted into
  increasing row order.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("dgBCMatrix",
    ...)}.

  Objects in this class have a structure similar to those in the
  \code{dgCMatrix} class except that each nonzero \dQuote{element} is
  itself a dense matrix with a given number of rows and columns.  These
  objects are used in the representation of linear mixed-effects model
  structures in the \code{lme4} package.
}
\section{Slots}{
  \describe{
    \item{\code{p}:}{Object of class \code{"integer"} of pointers, one
      for each column, to the initial (zero-based) index of elements in
      the column.}
    \item{\code{i}:}{Object of class \code{"integer"} of length nnzero
      (number of non-zero elements).  These are the row numbers for
      each non-zero element in the matrix.}
    \item{\code{x}:}{Object of class \code{"array"} - the non-zero
      elements of the matrix.  This is a three-dimensional array with
      the third dimension being nnzero.}
  }
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "dgBCMatrix", to = "dgCMatrix")}}
    \item{coerce}{\code{signature(from = "dgBCMatrix", to = "dgTMatrix")}}
    \item{coerce}{\code{signature(from = "dgBCMatrix", to = "dgeMatrix")}}
    \item{coerce}{\code{signature(from = "dgBCMatrix", to = "matrix")}}
    \item{coerce}{\code{signature(from = "dgCMatrix", to = "dgBCMatrix")}}
  }
}
%\references{}
%\author{}
%\note{}
\seealso{
 \code{\link{dgCMatrix-class}}
}
% \examples{

% }
\keyword{classes}
