\name{trMatrix-class}
\docType{class}
\alias{\%*\%,geMatrix,trMatrix-method}
\alias{\%*\%,trMatrix,geMatrix-method}
\alias{coerce,trMatrix,geMatrix-method}
\alias{coerce,trMatrix,matrix-method}
\alias{crossprod,trMatrix,missing-method}
\alias{determinant,trMatrix,logical-method}
\alias{determinant,trMatrix,missing-method}
\alias{norm,trMatrix,character-method}
\alias{norm,trMatrix,missing-method}
\alias{rcond,trMatrix,character-method}
\alias{rcond,trMatrix,missing-method}
\alias{solve,trMatrix,matrix-method}
\alias{solve,trMatrix,missing-method}
\alias{t,trMatrix-method}
\alias{trMatrix-class}
\title{Triangular, dense, non-packed, real matrices}
\description{The \code{"trMatrix"} class is the class of triangular
  matrices in dense, non-packed storage.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("trMatrix", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{uplo}:}{Object of class \code{"character"}. Must be
      either "U", for upper triangular, and "L", for lower triangular.}
    \item{\code{diag}:}{Object of class \code{"character"}. Must be
      either \code{"U"}, for unit triangular (diagonal is all ones), or
      \code{"N"} for non-unit.  The implicit diagonal elements are not
      explicitly stored when \code{diag} is \code{"U"}.}
    \item{\code{x}:}{Object of class \code{"numeric"}. The numeric
      values that constitute the matrix, stored in column-major order.}
    \item{\code{Dim}:}{Object of class \code{"integer"}. The dimensions
      of the matrix which must be a two-element vector of non-negative
      integers.}
    \item{\code{norm}:}{Object of class \code{"numeric"}. A named
      numeric vector of norms of the matrix.  Some possible components
      are \code{"O"}, the one norm; \code{"I"}, the infinity norm;
      \code{"F"}, the Frobenius norm; and \code{"M"}, the maximum
      absolute value.}
    \item{\code{rcond}:}{Object of class \code{"numeric"}. A named
      numeric vector of reciprocal condition numbers in either the
      1-norm \code{"O"} or the infinity norm \code{"I"}.}
    \item{\code{factorization}:}{Object of class \code{"list"}.  A named
      list of factorizations that have been computed for the matrix.}
  }
}
\section{Extends}{
Class \code{"geMatrix"}, directly.
Class \code{"Matrix"}, by class \code{"geMatrix"}.
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "trMatrix", to = "geMatrix")}}
    \item{coerce}{\code{signature(from = "trMatrix", to = "matrix")}}
    \item{norm}{\code{signature(x = "trMatrix", type = "character")}}
    \item{rcond}{\code{signature(x = "trMatrix", type = "character")}}    
    \item{solve}{\code{signature(a = "trMatrix", b = "missing")}}
    \item{solve}{\code{signature(a = "trMatrix", b = "matrix")}}
  }
}
%\references{}
%\author{}
\seealso{
  \code{\link{geMatrix-class}}, \code{\link{Matrix-class}}
}
%\examples{}
\keyword{classes}
