\name{dgeMatrix-class}
\title{General S4 Matrix class}
\docType{class}
\alias{dgeMatrix-class}
% FIXME: These two are not (yet) further referenced below:
\alias{dMatrix-class}
\alias{ddenseMatrix-class}
% METHODS:
% Group ones
\alias{Arith,dgeMatrix,dgeMatrix-method}
\alias{Arith,dgeMatrix,numeric-method}
\alias{Arith,numeric,dgeMatrix-method}
\alias{Math,dgeMatrix-method}
\alias{Math2,dgeMatrix,numeric-method}
% others
\alias{\%*\%,dgeMatrix,dgeMatrix-method}
\alias{Schur,dgeMatrix,logical-method}
\alias{Schur,dgeMatrix,missing-method}
\alias{coerce,matrix,dgeMatrix-method}
\alias{coerce,dgeMatrix,matrix-method}
\alias{crossprod,dgeMatrix,dgeMatrix-method}
\alias{crossprod,dgeMatrix,matrix-method}
\alias{crossprod,dgeMatrix,missing-method}
\alias{crossprod,dgeMatrix,numeric-method}
\alias{determinant,dgeMatrix,logical-method}
\alias{determinant,dgeMatrix,missing-method}
\alias{diag,dgeMatrix-method}
\alias{dim,dgeMatrix-method}
\alias{dimnames,dgeMatrix-method}
\alias{dimnames<-,dgeMatrix,list-method}
\alias{eigen,dgeMatrix,missing-method}
\alias{eigen,dgeMatrix,logical-method}
\alias{lu,dgeMatrix-method}
\alias{norm,dgeMatrix,missing-method}
\alias{norm,dgeMatrix,character-method}
\alias{rcond,dgeMatrix,missing-method}
\alias{rcond,dgeMatrix,character-method}
\alias{solve,dgeMatrix,matrix-method}
\alias{solve,dgeMatrix,missing-method}
\alias{t,dgeMatrix-method}
%
\description{A general, numeric matrix in the S4 Matrix representation}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("dgeMatrix", ...)}
  or, more commonly, by coercion from the \code{Matrix} class (see
  \link{Matrix-class}).
}
\section{Slots}{
  \describe{
    \item{\code{x}:}{Object of class \code{"numeric"} - the numeric
      values contained in the matrix, in column-major order.}
    \item{\code{Dim}:}{Object of class \code{"integer"} - the dimensions
     of the matrix - must be an integer vector with exactly two
     non-negative values.}
    \item{\code{Dimnames}:}{a list of length two - inherited from class
      \code{Matrix}, see \code{\link{Matrix-class}}.}
    \item{\code{rcond}:}{Object of class \code{"numeric"} - an estimate
      of the reciprocal of the condition number, if it has been computed
      at some point.}
    \item{\code{factors}:}{Object of class \code{"list"} - a list
      of factorizations of the matrix. }
  }
}
\section{Methods}{
  The are group methods (see, e.g., \code{\link{Arith}})
  \describe{
    \item{Arith}{\code{signature(e1 = "dgeMatrix", e2 = "dgeMatrix")}: ... }
    \item{Arith}{\code{signature(e1 = "dgeMatrix", e2 = "numeric")}: ... }
    \item{Arith}{\code{signature(e1 = "numeric", e2 = "dgeMatrix")}: ... }
    \item{Math}{\code{signature(x = "dgeMatrix")}: ... }
    \item{Math2}{\code{signature(x = "dgeMatrix", digits = "numeric")}: ...}
  }
  and other matrix methods available:
  \describe{
    \item{Schur}{\code{signature(x = "dgeMatrix", vectors = "logical")}: ... }
    \item{Schur}{\code{signature(x = "dgeMatrix", vectors = "missing")}: ... }
    \item{coerce}{\code{signature(from = "dgeMatrix", to = "matrix")}: ... }
    \item{coerce}{\code{signature(from = "matrix", to = "dgeMatrix")}: ... }
    \item{crossprod}{\code{signature(x = "dgeMatrix", y = "dgeMatrix")}: ... }
    \item{crossprod}{\code{signature(x = "dgeMatrix", y = "matrix")}: ... }
    \item{crossprod}{\code{signature(x = "dgeMatrix", y = "missing")}: ... }
    \item{crossprod}{\code{signature(x = "dgeMatrix", y = "numeric")}: ... }
    \item{diag}{\code{signature(x = "dgeMatrix")}: ... }
    \item{dim}{\code{signature(x = "dgeMatrix")}: ... }
    \item{dimnames}{\code{signature(x = "dgeMatrix")}: ... }
    \item{dimnames<-}{\code{signature(x = "dgeMatrix", value = "list")}: ... }
    \item{eigen}{\code{signature(x = "dgeMatrix", only.values= "logical")}: ...}
    \item{eigen}{\code{signature(x = "dgeMatrix", only.values= "missing")}: ...}
    \item{lu}{\code{signature(x = "dgeMatrix")}: ... }
    \item{norm}{\code{signature(x = "dgeMatrix", type = "character")}: ... }
    \item{norm}{\code{signature(x = "dgeMatrix", type = "missing")}: ... }
    \item{rcond}{\code{signature(x = "dgeMatrix", type = "character")}: ... }
    \item{rcond}{\code{signature(x = "dgeMatrix", type = "missing")}: ... }
    \item{solve}{\code{signature(a = "dgeMatrix", b = "dgeMatrix")}: ... }
    \item{solve}{\code{signature(a = "dgeMatrix", b = "missing")}: ... }
    \item{t}{\code{signature(x = "dgeMatrix")}: ... }
  }
}
%\references{}
%\author{}
%\note{}
\seealso{
  \code{\link{Matrix-class}},
  \code{\link{dtrMatrix-class}},
  \code{\link{dsyMatrix-class}}
}
%\examples{}
\keyword{classes}
\keyword{algebra}
