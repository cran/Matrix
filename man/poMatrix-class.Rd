\name{poMatrix-class}
\docType{class}
\alias{poMatrix-class}
\alias{rcond,poMatrix,character-method}
\alias{rcond,poMatrix,missing-method}
\alias{chol,poMatrix-method}
\alias{chol,poMatrix,ANY-method}
\alias{determinant,poMatrix,logical-method}
\alias{solve,poMatrix,geMatrix-method}
\alias{solve,poMatrix,matrix-method}
\alias{solve,poMatrix,missing-method}
\alias{solve,poMatrix,numeric-method}
\title{Positive semi-definite, dense, real matrices}
\description{The \code{"poMatrix"} class is the class of
  positive-semidefinite symmetric matrices in non-packed storage. Only the
  upper triangle or the lower triangle is required to be available.}
\section{Objects from the Class}{Objects can be created by calls of the
  form \code{new("poMatrix", ...)} or from \code{crossprod} applied to
  an \code{"geMatrix"} object.}
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
    \item{\code{factorization}:}{Object of class \code{"list"}.  A named
      list of factorizations that have been computed for the matrix.}
  }
}
\section{Extends}{
Class \code{"syMatrix"}, directly.
Class \code{"geMatrix"}, by class \code{"syMatrix"}.
Class \code{"Matrix"}, by class \code{"syMatrix"}.
}
\section{Methods}{
  \describe{
    \item{chol}{\code{signature(x = "poMatrix")}:
      Returns (and stores) the Cholesky decomposition of the matrix
      \code{x}.}
    \item{rcond}{\code{signature(x = "poMatrix", type = "character")}:
      Returns (and stores) the reciprocal of the condition number of
      \code{x}.  The \code{type} can be \code{"O"} for the
      one-norm (the default) or \code{"I"} for the infinity-norm.  For
      symmetric matrices the result does not depend on the type.}
    \item{solve}{\code{signature(a = "poMatrix", b = "missing")}:
      Return the inverse of \code{a}.}
    \item{solve}{\code{signature(a = "poMatrix", b = "numeric")}:
      Solve the linear system defined by \code{a} and \code{b}, where
      \code{b} can be a numeric vector, or a matrix, or a geMatrix
      object.  The Cholesky decomposition of \code{a} is calculated (if
      needed) while solving the system.}
  }
}
%\references{}
%\author{}
\seealso{
  \code{\link{syMatrix-class}}, \code{\link{geMatrix-class}},
  \code{\link{Matrix}}, \code{\link{rcond}}, \code{\link{chol}},
  \code{\link{solve}}, \code{\link{crossprod}}
}
%\examples{}
\keyword{classes}
\keyword{algebra}
