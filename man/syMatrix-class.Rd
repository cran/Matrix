\name{syMatrix-class}
\docType{class}
\alias{%*%,geMatrix,syMatrix-method}
\alias{%*%,syMatrix,geMatrix-method}
\alias{syMatrix-class}
\alias{coerce,syMatrix,geMatrix-method}
\alias{coerce,syMatrix,matrix-method}
\alias{norm,syMatrix,character-method}
\alias{norm,syMatrix,missing-method}
\alias{t,syMatrix-method}
\title{Class "syMatrix" of symmetric matrices}
\description{The \code{"syMatrix"} class is the class of symmetric
  matrices in non-packed storage.  Only the upper triangle or the lower
  triangle is stored.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("syMatrix", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{uplo}:}{Object of class \code{"character"}. Must be
      either "U", for upper triangular, and "L", for lower triangular.}
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
    \item{coerce}{\code{signature(from = "syMatrix", to = "geMatrix")}}
    \item{coerce}{\code{signature(from = "syMatrix", to = "matrix")}}
    \item{norm}{\code{signature(x = "syMatrix", type = "character")}}
  }
}
%\references{}
%\author{}
\seealso{
  \code{\link{geMatrix-class}}, \code{\link{Matrix-class}}
}
%\examples{}
\keyword{classes}
