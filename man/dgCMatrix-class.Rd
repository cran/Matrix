\name{dgCMatrix-class}
\docType{class}
\alias{dgCMatrix-class}
\alias{\%*\%,dgCMatrix,dgeMatrix-method}
\alias{\%*\%,dgCMatrix,matrix-method}
\alias{coerce,matrix,dgCMatrix-method}
\alias{coerce,dgCMatrix,dgTMatrix-method}
\alias{coerce,dgCMatrix,dgBCMatrix-method}
\alias{coerce,dgCMatrix,matrix-method}
\alias{coerce,dgCMatrix,dgeMatrix-method}
\alias{crossprod,dgCMatrix,missing-method}
\alias{crossprod,dgCMatrix,matrix-method}
\alias{crossprod,dgCMatrix,dgeMatrix-method}
\alias{diag,dgCMatrix-method}
\alias{dim,dgCMatrix-method}
\alias{image,dgCMatrix-method}
\alias{t,dgCMatrix-method}
\title{Compressed, sparse, column-oriented numeric matrices}
\description{The \code{dgCMatrix} class is a class of sparse numeric
  matrices in the compressed, sparse, column-oriented format.  In this
  implementation the non-zero elements in the columns are sorted into
  increasing row order.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("dgCMatrix", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{p}:}{Object of class \code{"integer"} of pointers, one
      for each column, to the initial (zero-based) index of elements in
      the column.}
    \item{\code{i}:}{Object of class \code{"integer"} of length nnzero
      (number of non-zero elements).  These are the row numbers for
      each non-zero element in the matrix.}
    \item{\code{x}:}{Object of class \code{"numeric"} - the non-zero
      elements of the matrix.}
    \item{\code{nrow}:}{Object of class \code{"integer"} - the number of
      rows in the matrix.}
  }
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "matrix", to = "dgCMatrix")}}
    \item{coerce}{\code{signature(from = "dgCMatrix", to = "matrix")}}
    \item{coerce}{\code{signature(from = "dgCMatrix", to = "dgBCMatrix")}}
    \item{coerce}{\code{signature(from = "dgCMatrix", to = "dgTMatrix")}}
    \item{crossprod}{\code{signature(x = "dgCMatrix", y = "missing")}
      returns \code{t(x) \%*\% x} as an \code{dsCMatrix} object.}
    \item{crossprod}{\code{signature(x = "dgCMatrix", y = "matrix")}: ... }
    \item{crossprod}{\code{signature(x = "dgCMatrix", y = "numeric")}: ... }
    \item{diag}{\code{signature(x = "dgCMatrix")}: returns the diagonal
      of \code{x}}
    \item{dim}{\code{signature(x = "dgCMatrix")}: returns the dimensions
      of \code{x}}
    \item{image}{\code{signature(x = "dgCMatrix")}: plots an image of
      \code{x} using the \code{\link[lattice]{levelplot}} function}
  }
}
%\references{}
%\author{}
%\note{}
\seealso{
  \code{\link{dsCMatrix-class}}, \code{\link{dtCMatrix-class}}
}
%\examples{}
\keyword{classes}
\keyword{algebra}
