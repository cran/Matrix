\name{dsyMatrix-class}
\docType{class}
\alias{\%*\%,dgeMatrix,dsyMatrix-method}
\alias{\%*\%,dsyMatrix,dgeMatrix-method}
\alias{dsyMatrix-class}
\alias{dspMatrix-class}
\alias{coerce,dsyMatrix,dgeMatrix-method}
\alias{coerce,dsyMatrix,matrix-method}
\alias{norm,dsyMatrix,character-method}
\alias{norm,dsyMatrix,missing-method}
\alias{t,dsyMatrix-method}
\title{Symmetric, dense matrices}
\description{The \code{"dsyMatrix"} class is the class of symmetric, dense
  matrices in non-packed storage and \code{"dsyMatrix"} is the class of
  symmetric dense matrices in packed storage.  Only the upper triangle or the lower
  triangle is stored.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("dsyMatrix", ...)}.
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
    \item{\code{rcond}:}{Object of class \code{"numeric"}. A named
      numeric vector of reciprocal condition numbers in either the
      1-norm \code{"O"} or the infinity norm \code{"I"}.}
    \item{\code{factors}:}{Object of class \code{"list"}.  A named
      list of factorizations that have been computed for the matrix.}
  }
}
\section{Extends}{
Class \code{"dgeMatrix"}, directly.
Class \code{"Matrix"}, by class \code{"dgeMatrix"}.
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "dsyMatrix", to = "dgeMatrix")}}
    \item{coerce}{\code{signature(from = "dsyMatrix", to = "matrix")}}
    \item{norm}{\code{signature(x = "dsyMatrix", type = "character")}}
  }
}
%\references{}
%\author{}
\seealso{
  \code{\link{dgeMatrix-class}}, \code{\link{Matrix-class}}
}
%\examples{}
\keyword{classes}
