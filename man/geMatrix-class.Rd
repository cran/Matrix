\name{geMatrix-class}
\docType{class}
\alias{geMatrix-class}
\alias{%*%,geMatrix,geMatrix-method}
\alias{Schur,geMatrix,logical-method}
\alias{Schur,geMatrix,missing-method}
\alias{coerce,matrix,geMatrix-method}
\alias{coerce,geMatrix,matrix-method}
\alias{crossprod,geMatrix,geMatrix-method}
\alias{crossprod,geMatrix,matrix-method}
\alias{crossprod,geMatrix,missing-method}
\alias{crossprod,geMatrix,numeric-method}
\alias{determinant,geMatrix,logical-method}
\alias{determinant,geMatrix,missing-method}
\alias{diag,geMatrix-method}
\alias{dim,geMatrix-method}
\alias{lu,geMatrix-method}
\alias{norm,geMatrix,missing-method}
\alias{norm,geMatrix,character-method}
\alias{rcond,geMatrix,missing-method}
\alias{rcond,geMatrix,character-method}
\alias{solve,geMatrix,matrix-method}
\alias{solve,geMatrix,missing-method}
\title{General S4 Matrix class}
\description{A general, numeric matrix in the S4 Matrix representation}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("geMatrix", ...)}
or, more commonly, by coercion from the \code{matrix} class.
}
\section{Slots}{
  \describe{
    \item{\code{x}:}{Object of class \code{"numeric"} - the numeric
      values contained in the matrix, in column-major order.}
    \item{\code{Dim}:}{Object of class \code{"integer"} - the dimensions
     of the matrix - must be an integer vector with exactly two
     non-negative values.}
    \item{\code{rcond}:}{Object of class \code{"numeric"} - an estimate
      of the reciprocal of the condition number, if it has been computed
      at some point.}
    \item{\code{factorization}:}{Object of class \code{"list"} - a list
      of factorizations of the matrix. }
  }
}
\section{Methods}{
  \describe{
    \item{Schur}{\code{signature(x = "geMatrix", vectors = "logical")}: ... }
    \item{Schur}{\code{signature(x = "geMatrix", vectors = "missing")}: ... }
    \item{coerce}{\code{signature(from = "geMatrix", to = "matrix")}: ... }
    \item{coerce}{\code{signature(from = "matrix", to = "geMatrix")}: ... }
    \item{coerce}{\code{signature(from = "matrix", to = "geMatrix")}}
    \item{crossprod}{\code{signature(x = "geMatrix", y = "geMatrix")}: ... }
    \item{crossprod}{\code{signature(x = "geMatrix", y = "matrix")}: ... }
    \item{crossprod}{\code{signature(x = "geMatrix", y = "missing")}: ... }
    \item{crossprod}{\code{signature(x = "geMatrix", y = "numeric")}: ... }
    \item{diag}{\code{signature(x = "geMatrix")}: ... }
    \item{dim}{\code{signature(x = "geMatrix")}: ... }
    \item{lu}{\code{signature(x = "geMatrix")}: ... }
    \item{norm}{\code{signature(x = "geMatrix", type = "character")}: ... }
    \item{norm}{\code{signature(x = "geMatrix", type = "missing")}: ... }
    \item{rcond}{\code{signature(x = "geMatrix", type = "character")}: ... }
    \item{rcond}{\code{signature(x = "geMatrix", type = "missing")}: ... }
    \item{solve}{\code{signature(a = "geMatrix", b = "geMatrix")}: ... }
    \item{solve}{\code{signature(a = "geMatrix", b = "missing")}: ... }
  }
}
%\references{}
%\author{}
%\note{}
\seealso{
  \code{\link{Matrix-class}}, \code{\link{trMatrix-class}}, \code{\link{syMatrix-class}}
}
%\examples{}
\keyword{classes}
\keyword{algebra}
