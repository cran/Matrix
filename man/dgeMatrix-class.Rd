\name{dgeMatrix-class}
\title{Class "dgeMatrix" of Dense Numeric (S4 Class) Matrices}
\docType{class}
%
\alias{dgeMatrix-class}
%
\alias{Arith,dgeMatrix,dgeMatrix-method}
\alias{Arith,dgeMatrix,logical-method}
\alias{Arith,dgeMatrix,numeric-method}
\alias{Arith,dgeMatrix,sparseVector-method}
\alias{Arith,logical,dgeMatrix-method}
\alias{Arith,numeric,dgeMatrix-method}
\alias{Math,dgeMatrix-method}
\alias{as.vector,dgeMatrix-method}
\alias{coerce,dgeMatrix,matrix-method}
\alias{coerce,dgeMatrix,vector-method}
\alias{determinant,dgeMatrix,logical-method}
\alias{log,dgeMatrix-method}
%
\description{A general numeric dense matrix in the S4 Matrix
  representation.  \code{dgeMatrix} is the \emph{\dQuote{standard}}
  class for dense numeric matrices in the \pkg{Matrix} package.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("dgeMatrix", ...)}
  or, more commonly, by coercion from the \code{Matrix} class (see
  \linkS4class{Matrix}) or by \code{\link{Matrix}(..)}.
}
\section{Slots}{
  \describe{
    \item{\code{x}:}{Object of class \code{"numeric"} - the numeric
      values contained in the matrix, in column-major order.}
    \item{\code{Dim}:}{Object of class \code{"integer"} - the dimensions
     of the matrix - must be an integer vector with exactly two
     non-negative values.}
    \item{\code{Dimnames}:}{a list of length two - inherited from class
      \code{\linkS4class{Matrix}}.}
    \item{\code{factors}:}{Object of class \code{"list"} - a list
      of factorizations of the matrix.}
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
  matrix products \code{\link[=crossprod-methods]{\%*\%}},
  \code{\link[=crossprod-methods]{crossprod}()} and \code{tcrossprod()},
  several \code{\link{solve}} methods,
  and other matrix methods available:
  \describe{
    \item{Schur}{\code{signature(x = "dgeMatrix", vectors = "logical")}: ... }
    \item{Schur}{\code{signature(x = "dgeMatrix", vectors = "missing")}: ... }
    \item{chol}{\code{signature(x = "dgeMatrix")}: see \code{\link{chol}}.}
    \item{colMeans}{\code{signature(x = "dgeMatrix")}: columnwise means (averages)}
    \item{colSums}{\code{signature(x = "dgeMatrix")}: columnwise sums}
    \item{diag}{\code{signature(x = "dgeMatrix")}: ... }
    \item{dim}{\code{signature(x = "dgeMatrix")}: ... }
    \item{dimnames}{\code{signature(x = "dgeMatrix")}: ... }
    \item{eigen}{\code{signature(x = "dgeMatrix", only.values= "logical")}: ...}
    \item{eigen}{\code{signature(x = "dgeMatrix", only.values= "missing")}: ...}
    \item{norm}{\code{signature(x = "dgeMatrix", type = "character")}: ... }
    \item{norm}{\code{signature(x = "dgeMatrix", type = "missing")}: ... }
    \item{rcond}{\code{signature(x = "dgeMatrix", norm = "character")}
      or \code{norm = "missing"}:
      the reciprocal condition number, \code{\link{rcond}()}.}
    \item{rowMeans}{\code{signature(x = "dgeMatrix")}: rowwise means (averages)}
    \item{rowSums}{\code{signature(x = "dgeMatrix")}: rowwise sums}
    \item{t}{\code{signature(x = "dgeMatrix")}: matrix transpose}
  }
}
\seealso{
  Classes \code{\linkS4class{Matrix}},
  \code{\linkS4class{dtrMatrix}}, and \code{\linkS4class{dsyMatrix}}.

}
%\examples{}
\keyword{classes}
\keyword{algebra}
