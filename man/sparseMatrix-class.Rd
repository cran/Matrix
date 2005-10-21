\name{sparseMatrix-class}
\docType{class}
\title{Virtual Class "sparseMatrix" --- Mother of Sparse Matrices}
\alias{sparseMatrix-class}

\alias{-,sparseMatrix,missing-method}
\alias{show,sparseMatrix-method}
\alias{coerce,graph,sparseMatrix-method}
\alias{coerce,graphNEL,sparseMatrix-method}
\alias{coerce,sparseMatrix,graph-method}
\alias{coerce,sparseMatrix,graphNEL-method}
\alias{coerce,dgTMatrix,graphNEL-method}
\description{Virtual Mother Class of All Sparse Matrices}
\section{Slots}{
  \describe{
    \item{\code{Dim}:}{Object of class \code{"integer"} - the dimensions
     of the matrix - must be an integer vector with exactly two
     non-negative values.}
    \item{\code{Dimnames}:}{a list of length two - inherited from class
      \code{Matrix}, see \code{\link{Matrix-class}}.}
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
      \emph{\dQuote{structural}} zeroes as \code{"."} or blank.}
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
}
\keyword{classes}
