\name{sscMatrix-class}
\docType{class}
\alias{sscMatrix-class}
\alias{solve,sscMatrix,matrix-method}
\alias{chol,sscMatrix,missing-method}
\alias{chol,sscMatrix,logical-method}
\alias{t,sscMatrix-method}
\title{Symmetric, compressed, sparse column matrices}
\description{The \code{sscMatrix} class is a class of symmetric, sparse
  matrices in the compressed, column-oriented format.  In this
  implementation the non-zero elements in the columns are sorted into
  increasing row order.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("sscMatrix", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{uplo}:}{A character object indicating if the upper
      triangle (\code{"U"} or \code{"u"}) or the lower triangle
      (\code{"L"} or \code{"l"}) is stored.  At present only the lower
      triangle form is allowed.}
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
\section{Methods}{
  \describe{
    \item{solve}{\code{signature(a = "sscMatrix", b = "matrix")}: Solve
      a linear system of equations defined by \code{x} using a Cholesky
      decomposition.}
    \item{chol}{\code{signature(x = "sscMatrix", pivot = "logical")}:
      Returns (and stores) the Cholesky decomposition of the matrix
      \code{x}. If \code{pivot} is \code{TRUE} (the default) Metis is
      used to create a reordering of the rows and columns of \code{x} so
      as to minimize fill-in.}
  }
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
