\name{cscMatrix-class}
\docType{class}
\alias{cscMatrix-class}
\alias{coerce,matrix,cscMatrix-method}
\alias{coerce,cscMatrix,tripletMatrix-method}
\alias{coerce,cscMatrix,matrix-method}
\alias{coerce,cscMatrix,geMatrix-method}
\alias{crossprod,cscMatrix,missing-method}
\alias{crossprod,cscMatrix,matrix-method}
\alias{crossprod,cscMatrix,numeric-method}
\alias{diag,cscMatrix-method}
\alias{dim,cscMatrix-method}
\alias{image,cscMatrix-method}
\alias{t,cscMatrix-method}
\title{Class "cscMatrix" compressed, sparse, column Matrix}
\description{The \code{cscMatrix} class is a class of sparse matrices in
  the compressed, sparse, column-oriented format.  In this
  implementation the non-zero elements in the columns are sorted into
  increasing row order.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("cscMatrix", ...)}.
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
    \item{coerce}{\code{signature(from = "matrix", to = "cscMatrix")}}
    \item{coerce}{\code{signature(from = "cscMatrix", to = "matrix")}}
    \item{coerce}{\code{signature(from = "cscMatrix", to = "tripletMatrix")}}
    \item{crossprod}{\code{signature(x = "cscMatrix", y = "missing")}
      returns \code{t(x) \%*\% x} as an \code{sscMatrix} object.}
    \item{crossprod}{\code{signature(x = "cscMatrix", y = "matrix")}: ... }
    \item{crossprod}{\code{signature(x = "cscMatrix", y = "numeric")}: ... }
    \item{diag}{\code{signature(x = "cscMatrix")}: returns the diagonal
      of \code{x}}
    \item{dim}{\code{signature(x = "cscMatrix")}: returns the dimensions
      of \code{x}}
    \item{image}{\code{signature(x = "cscMatrix")}: plots an image of
      \code{x} using the \code{\link[lattice]{levelplot}} function}
  }
}
%\references{}
%\author{}
%\note{}
\seealso{
  \code{\link{sscMatrix-class}}, \code{\link{tscMatrix-class}}
}
%\examples{}
\keyword{classes}
\keyword{algebra}
