\name{dtpMatrix-class}
\docType{class}
\alias{dtpMatrix-class}
\alias{\%*\%,dtpMatrix,dgeMatrix-method}
\alias{\%*\%,dgeMatrix,dtpMatrix-method}
\alias{\%*\%,dtpMatrix,dtpMatrix-method}
\alias{\%*\%,dtpMatrix,matrix-method}
\alias{\%*\%,matrix,dtpMatrix-method}
\alias{\%*\%,dtpMatrix,numeric-method}
\alias{\%*\%,numeric,dtpMatrix-method}
\alias{coerce,dtpMatrix,dtrMatrix-method}
\alias{coerce,dtpMatrix,dgeMatrix-method}
\alias{coerce,dtpMatrix,matrix-method}
\alias{coerce,matrix,dtpMatrix-method}
\alias{determinant,dtpMatrix,missing-method}
\alias{determinant,dtpMatrix,logical-method}
\alias{diag,dtpMatrix-method}
\alias{norm,dtpMatrix,character-method}
\alias{norm,dtpMatrix,missing-method}
\alias{rcond,dtpMatrix,character-method}
\alias{rcond,dtpMatrix,missing-method}
\alias{solve,dtpMatrix,missing-method}
\alias{solve,dtpMatrix,matrix-method}
\alias{t,dtpMatrix-method}
\alias{unpack,dtpMatrix-method}
\title{Packed triangular dense matrices}
\description{The \code{"dtpMatrix"} class is the class of triangular,
  dense, numeric matrices in packed storage.  The \code{"dtrMatrix"}
  class is the same except in nonpacked storage.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("dtpMatrix",
    ...)} or by coercion from other classes of matrices.
}
\section{Slots}{
  \describe{
    \item{\code{uplo}:}{Object of class \code{"character"}. Must be
      either "U", for upper triangular, and "L", for lower triangular.}
    \item{\code{diag}:}{Object of class \code{"character"}. Must be
      either \code{"U"}, for unit triangular (diagonal is all ones), or
      \code{"N"} for non-unit.  The implicit diagonal elements are not
      explicitly stored when \code{diag} is \code{"U"}.}
%%% FIXME: The above is not true -- they *are* stored (but values not considered)
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
    \item{\code{Dimnames}:}{Object of class \code{"list"}. The dimnames
      of the matrix which must be a list of two elements.  The elements
      should be either \code{NULL} or a character vector whose length
      matches the corresponding dimension.}
  }
}
\section{Extends}{
Class \code{"ddenseMatrix"}, directly.
Class \code{"dMatrix"}, by class \code{"ddenseMatrix"}.
Class \code{"Matrix"}, by class \code{"ddenseMatrix"}.
}
\section{Methods}{
  \describe{
    \item{\%*\%}{\code{signature(x = "dtpMatrix", y = "dgeMatrix")}:
      Matrix multiplication; dito for several other signature
      combinations, see \code{showMethods("\%*\%", class = "dtpMatrix")}.}
    \item{coerce}{\code{signature(from = "dtpMatrix", to = "dtrMatrix")}}
    \item{coerce}{\code{signature(from = "dtpMatrix", to = "dgeMatrix")}}
    \item{coerce}{\code{signature(from = "dtpMatrix", to = "matrix")}}
    \item{determinant}{\code{signature(x = "dtpMatrix", logarithm = "missing")}: ... }
    \item{determinant}{\code{signature(x = "dtpMatrix", logarithm = "logical")}: ... }
    \item{diag}{\code{signature(x = "dtpMatrix")}: ... }
    \item{norm}{\code{signature(x = "dtpMatrix", type = "character")}: ... }
    \item{norm}{\code{signature(x = "dtpMatrix", type = "missing")}: ... }
    \item{rcond}{\code{signature(x = "dtpMatrix", type = "character")}: ... }
    \item{rcond}{\code{signature(x = "dtpMatrix", type = "missing")}: ... }
    \item{solve}{\code{signature(a = "dtpMatrix", b = "missing")}: ... }
    \item{solve}{\code{signature(a = "dtpMatrix", b = "matrix")}: ... }
    \item{t}{\code{signature(x = "dtpMatrix")}: ... }
    \item{unpack}{\code{signature(x = "dtpMatrix")}: ... }
  }
}
\seealso{
  \code{\link{dtrMatrix-class}}
}
\examples{
example("dtrMatrix-class")
(p1 <- as(T3, "dtpMatrix"))
str(p1)
(pp <- as(T, "dtpMatrix"))
}
\keyword{classes}
