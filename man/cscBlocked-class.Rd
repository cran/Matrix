\name{cscBlocked-class}
\docType{class}
\alias{cscBlocked-class}
\alias{coerce,cscBlocked,cscMatrix-method}
\alias{coerce,cscBlocked,tripletMatrix-method}
\title{Class "cscBlocked" compressed, sparse, blocked Matrix }
\description{The \code{cscMatrix} class is a class of sparse matrices in
  the compressed, sparse, column-oriented format.  In this
  implementation the non-zero elements in the columns are sorted into
  increasing row order.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("cscBlocked",
    ...)}.

  Objects in this class have a structure similar to those in the
  \code{cscMatrix} class except that each nonzero "element" is itself a
  dense matrix with a given number of rows and columns.  These objects
  are used in the representation of linear mixed-effects model
  structures, \code{lmer}.
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
    \item{coerce}{\code{signature(from = "cscBlocked", to = "cscMatrix")}}
    \item{coerce}{\code{signature(from = "cscBlocked", to = "tripletMatrix")}}
    \item{coerce}{\code{signature(from = "cscMatrix", to = "cscBlocked")}}
  }
}
%\references{}
%\author{}
%\note{}
\seealso{
 \code{\link{cscMatrix-class}}, \code{\link{lmer-class}}
}
%\examples{}
\keyword{classes}
