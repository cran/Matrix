\name{dtrMatrix-class}
\docType{class}
\alias{\%*\%,dtrMatrix,dtrMatrix-method}
\alias{\%*\%,dgeMatrix,dtrMatrix-method}
\alias{\%*\%,dtrMatrix,dgeMatrix-method}
\alias{\%*\%,dtrMatrix,matrix-method}
\alias{\%*\%,matrix,dtrMatrix-method}
\alias{coerce,dtrMatrix,dgeMatrix-method}
\alias{coerce,dtrMatrix,matrix-method}
\alias{coerce,dtrMatrix,dtpMatrix-method}
\alias{coerce,matrix,dtrMatrix-method}
\alias{crossprod,dtrMatrix,missing-method}
\alias{determinant,dtrMatrix,logical-method}
\alias{determinant,dtrMatrix,missing-method}
\alias{norm,dtrMatrix,character-method}
\alias{norm,dtrMatrix,missing-method}
\alias{rcond,dtrMatrix,character-method}
\alias{rcond,dtrMatrix,missing-method}
\alias{solve,dtrMatrix,dgeMatrix-method}
\alias{solve,dtrMatrix,matrix-method}
\alias{solve,dtrMatrix,missing-method}
\alias{t,dtrMatrix-method}
\alias{dtrMatrix-class}
\title{Triangular, dense, numeric matrices}
\description{
  The \code{"dtrMatrix"} class is the class of triangular, dense,
  numeric matrices in nonpacked storage.  The \code{"dtpMatrix"} class
  is the same except in packed storage.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("dtrMatrix", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{uplo}:}{Object of class \code{"character"}. Must be
      either "U", for upper triangular, and "L", for lower triangular.}
    \item{\code{diag}:}{Object of class \code{"character"}. Must be
      either \code{"U"}, for unit triangular (diagonal is all ones), or
      \code{"N"}; see \code{\link{triangularMatrix-class}}.}
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
Class \code{"triangularMatrix"}, directly.
Class \code{"Matrix"} and others, by class \code{"dgeMatrix"}.
}
\section{Methods}{
  \describe{
    \item{\%*\%}{\code{signature(x = "dtrMatrix", y = "matrix")} and other
      signatures (use \code{showMethods("\%*\%", class="dtrMatrix")}):
      matrix multiplication.}
    \item{coerce}{\code{signature(from = "dtrMatrix", to = "dgeMatrix")}}
    \item{coerce}{\code{signature(from = "dtrMatrix", to = "dgeMatrix")}}
    \item{coerce}{\code{signature(from = "dtrMatrix", to = "matrix")}}
    \item{coerce}{\code{signature(from = "dtrMatrix", to = "dtpMatrix")}}
    \item{norm}{\code{signature(x = "dtrMatrix", type = "character")}}
    \item{rcond}{\code{signature(x = "dtrMatrix", type = "character")}}
    \item{solve}{\code{signature(a = "dtrMatrix", b = "missing")}}
    \item{solve}{\code{signature(a = "dtrMatrix", b = "matrix")}}
  }
}
%\references{}
%\author{}
\seealso{
  \code{\link{dgeMatrix-class}}, \code{\link{dtpMatrix-class}},
  \code{\link{Matrix-class}}
}
\examples{
(m <- rbind(2:3, 0:-1))
(M <- as(m, "dgeMatrix"))
(T <- as(M, "dtrMatrix")) ## upper triangular is default
## hence this may surprise at first:
(T2 <- as(t(M), "dtrMatrix"))# is factually diagonal
## and different from
(T3 <- t(T)) # which is *lower* triangular:
str(T3)
}
\keyword{classes}
