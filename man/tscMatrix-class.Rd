\name{tscMatrix-class}
\docType{class}
\alias{tscMatrix-class}
\alias{t,tscMatrix-method}
\title{Triangular, compressed, sparse column matrices}
\description{The \code{sscMatrix} class is a class of triangular, sparse
  matrices in the compressed, column-oriented format.  In this
  implementation the non-zero elements in the columns are sorted into
  increasing row order.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("tscMatrix", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{uplo}:}{A character object indicating if the matrix is upper
      triangular (\code{"U"} or \code{"u"}) or lower triangular
      (\code{"L"} or \code{"l"}).}
    \item{\code{p}:}{Object of class \code{"integer"} of pointers, one
      for each column, to the initial (zero-based) index of elements in
      the column.}
    \item{\code{i}:}{Object of class \code{"integer"} of length nnzero
      (number of non-zero elements).  These are the row numbers for
      each non-zero element in the matrix.}
    \item{\code{x}:}{Object of class \code{"numeric"} - the non-zero
      elements of the matrix.}
    \item{\code{factorization}:}{Object of class \code{"list"} - a list
      of factorizations of the matrix. }
    \item{\code{Dim}:}{Object of class \code{"integer"} - the dimensions
     of the matrix - must be an integer vector with exactly two
     non-negative values.}
  }
}
\section{Extends}{
Class \code{"cscMatrix"}, directly.
}
%\references{}
%\author{}
%\note{}
\seealso{
  \code{\link{cscMatrix-class}}
}
%\examples{}
\keyword{classes}
\keyword{algebra}
