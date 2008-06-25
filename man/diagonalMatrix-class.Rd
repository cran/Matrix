\name{diagonalMatrix-class}
\title{Class "diagonalMatrix" of Diagonal Matrices}
\docType{class}
\alias{diagonalMatrix-class}
%
\alias{as.vector,diagonalMatrix,missing-method}
\alias{coerce,matrix,diagonalMatrix-method}
\alias{coerce,diagonalMatrix,generalMatrix-method}
\alias{coerce,diagonalMatrix,matrix-method}
\alias{coerce,diagonalMatrix,nMatrix-method}
\alias{coerce,Matrix,diagonalMatrix-method}
\alias{cbind2,sparseMatrix,diagonalMatrix-method}
\alias{cbind2,diagonalMatrix,sparseMatrix-method}
\alias{rbind2,sparseMatrix,diagonalMatrix-method}
\alias{rbind2,diagonalMatrix,sparseMatrix-method}
\alias{determinant,diagonalMatrix,logical-method}
\alias{norm,diagonalMatrix,character-method}
%
\alias{coerce,ddiMatrix,dgeMatrix-method}
\alias{coerce,ddiMatrix,ddenseMatrix-method}
\alias{coerce,ldiMatrix,ldenseMatrix-method}
% Mainly against ambiguity warnings:
% Horrible-Hack: currently define for "all subclasses of diagonalMatrix"
% ------------- in ../R/diagMatrix.R
\alias{coerce,ddiMatrix,symmetricMatrix-method}
\alias{coerce,ldiMatrix,symmetricMatrix-method}
\alias{coerce,ddiMatrix,triangularMatrix-method}
\alias{coerce,ldiMatrix,triangularMatrix-method}
\alias{coerce,ddiMatrix,sparseMatrix-method}
\alias{coerce,ldiMatrix,sparseMatrix-method}
\alias{coerce,ddiMatrix,CsparseMatrix-method}
\alias{coerce,ldiMatrix,CsparseMatrix-method}
\alias{coerce,ddiMatrix,TsparseMatrix-method}
\alias{coerce,ldiMatrix,TsparseMatrix-method}
%
\alias{cbind2,ddiMatrix,matrix-method}
\alias{cbind2,ldiMatrix,matrix-method}
\alias{cbind2,matrix,ddiMatrix-method}
\alias{cbind2,matrix,ldiMatrix-method}
\alias{rbind2,ddiMatrix,matrix-method}
\alias{rbind2,ldiMatrix,matrix-method}
\alias{rbind2,matrix,ddiMatrix-method}
\alias{rbind2,matrix,ldiMatrix-method}
%
\alias{cbind2,ddiMatrix,atomicVector-method}
\alias{cbind2,ldiMatrix,atomicVector-method}
\alias{cbind2,atomicVector,ddiMatrix-method}
\alias{cbind2,atomicVector,ldiMatrix-method}
\alias{rbind2,ddiMatrix,atomicVector-method}
\alias{rbind2,ldiMatrix,atomicVector-method}
\alias{rbind2,atomicVector,ddiMatrix-method}
\alias{rbind2,atomicVector,ldiMatrix-method}
\alias{diag,diagonalMatrix-method}
% \alias{diag,ddiMatrix-method}
% \alias{diag,ldiMatrix-method}
%
\alias{Ops,ddiMatrix,sparseMatrix-method}
\alias{Ops,sparseMatrix,ddiMatrix-method}
\alias{Ops,ldiMatrix,sparseMatrix-method}
\alias{Ops,sparseMatrix,ldiMatrix-method}
\alias{Ops,ddiMatrix,numeric-method}
\alias{Ops,numeric,ddiMatrix-method}
\alias{Ops,ldiMatrix,numeric-method}
\alias{Ops,numeric,ldiMatrix-method}
%
\alias{Ops,ddiMatrix,ANY-method}
\alias{Ops,ANY,ddiMatrix-method}
\alias{Ops,ldiMatrix,ANY-method}
\alias{Ops,ANY,ldiMatrix-method}
\alias{Ops,ddiMatrix,Matrix-method}
\alias{Ops,Matrix,ddiMatrix-method}
\alias{Ops,ldiMatrix,Matrix-method}
\alias{Ops,Matrix,ldiMatrix-method}
\alias{Ops,ddiMatrix,dMatrix-method}
\alias{Ops,dMatrix,ddiMatrix-method}
\alias{Ops,ldiMatrix,dMatrix-method}
\alias{Ops,dMatrix,ldiMatrix-method}
%
\alias{Ops,ddiMatrix,ddiMatrix-method}
\alias{Ops,ddiMatrix,ldiMatrix-method}
\alias{Ops,ldiMatrix,ddiMatrix-method}
\alias{Ops,ldiMatrix,ldiMatrix-method}
\alias{Arith,numeric,ddiMatrix-method}
\alias{Arith,numeric,ldiMatrix-method}
\alias{Arith,ddiMatrix,numeric-method}
\alias{Arith,ldiMatrix,numeric-method}

\alias{all,ddiMatrix-method}
\alias{all,ldiMatrix-method}
\alias{any,ddiMatrix-method}
\alias{any,ldiMatrix-method}
\alias{prod,ddiMatrix-method}
\alias{prod,ldiMatrix-method}
\alias{sum,ddiMatrix-method}
\alias{sum,ldiMatrix-method}

\alias{&,ANY,ddiMatrix-method}
\alias{&,ANY,ldiMatrix-method}
\alias{&,ddiMatrix,ANY-method}
\alias{&,ldiMatrix,ANY-method}

\alias{&,Matrix,ddiMatrix-method}
\alias{&,Matrix,ldiMatrix-method}
\alias{&,ddiMatrix,Matrix-method}
\alias{&,ldiMatrix,Matrix-method}

\alias{&,ddiMatrix,denseMatrix-method}
\alias{&,ldiMatrix,denseMatrix-method}
\alias{&,denseMatrix,ddiMatrix-method}
\alias{&,denseMatrix,ldiMatrix-method}
%-end{H.Hack}
\alias{Ops,diagonalMatrix,diagonalMatrix-method}
\alias{Ops,diagonalMatrix,ddiMatrix-method}
\alias{Ops,ddiMatrix,diagonalMatrix-method}

\alias{Summary,ddiMatrix-method}
\alias{Summary,ldiMatrix-method}
%
\alias{t,diagonalMatrix-method}
\alias{show,diagonalMatrix-method}
\alias{\%*\%,dgeMatrix,diagonalMatrix-method}
\alias{\%*\%,matrix,diagonalMatrix-method}
\alias{\%*\%,diagonalMatrix,dgeMatrix-method}
\alias{\%*\%,diagonalMatrix,diagonalMatrix-method}
\alias{\%*\%,diagonalMatrix,matrix-method}
\alias{\%*\%,diagonalMatrix,CsparseMatrix-method}
\alias{\%*\%,diagonalMatrix,sparseMatrix-method}
\alias{\%*\%,CsparseMatrix,diagonalMatrix-method}
\alias{\%*\%,sparseMatrix,diagonalMatrix-method}
\alias{\%*\%,denseMatrix,diagonalMatrix-method}
\alias{\%*\%,diagonalMatrix,denseMatrix-method}
\alias{\%*\%,diagonalMatrix,lgeMatrix-method}
\alias{\%*\%,lgeMatrix,diagonalMatrix-method}

\alias{crossprod,diagonalMatrix,dgeMatrix-method}
\alias{crossprod,diagonalMatrix,lgeMatrix-method}
\alias{crossprod,diagonalMatrix,diagonalMatrix-method}
\alias{crossprod,diagonalMatrix,matrix-method}
\alias{crossprod,diagonalMatrix,missing-method}
\alias{crossprod,diagonalMatrix,sparseMatrix-method}
\alias{crossprod,sparseMatrix,diagonalMatrix-method}
\alias{tcrossprod,diagonalMatrix,diagonalMatrix-method}
\alias{tcrossprod,diagonalMatrix,sparseMatrix-method}
\alias{tcrossprod,diagonalMatrix,matrix-method}
\alias{tcrossprod,matrix,diagonalMatrix-method}
\alias{tcrossprod,dgeMatrix,diagonalMatrix-method}
\alias{tcrossprod,lgeMatrix,diagonalMatrix-method}
\alias{tcrossprod,diagonalMatrix,missing-method}
\alias{tcrossprod,sparseMatrix,diagonalMatrix-method}
\alias{solve,diagonalMatrix,Matrix-method}
\alias{solve,diagonalMatrix,matrix-method}
\alias{solve,diagonalMatrix,missing-method}
\description{
  Class "diagonalMatrix" is the virtual class of all diagonal matrices.
}
\section{Objects from the Class}{A virtual Class: No objects may be
  created from it.}
\section{Slots}{
  \describe{
    \item{\code{diag}:}{code{"character"} string, either \code{"U"} or
      \code{"N"}, where \code{"U"} means \sQuote{unit-diagonal}.}
    \item{\code{Dim}:}{matrix dimension, and}
    \item{\code{Dimnames}:}{the \code{\link{dimnames}}, a
      \code{\link{list}}, see the \code{\linkS4class{Matrix}} class
      description.  Typically \code{list(NULL,NULL)} for diagonal matrices.}
  }
}
\section{Extends}{
  Class \code{"\linkS4class{denseMatrix}"}, directly.
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "matrix", to = "diagonalMatrix")}: ... }
    \item{coerce}{\code{signature(from = "diagonalMatrix", to = "symmetricMatrix")}: ... }
    \item{coerce}{\code{signature(from = "diagonalMatrix", to = "triangularMatrix")}: ... }
    \item{coerce}{\code{signature(from = "diagonalMatrix", to = "matrix")}: ... }
    \item{coerce}{\code{signature(from = "Matrix", to = "diagonalMatrix")}: ... }
    \item{t}{\code{signature(x = "diagonalMatrix")}: ... }
    \item{\%*\%}{\code{signature(x = "dgeMatrix", y = "diagonalMatrix")}: ... }
    \item{\%*\%}{\code{signature(x = "matrix", y = "diagonalMatrix")}: ... }
    \item{\%*\%}{\code{signature(x = "diagonalMatrix", y = "matrix")}: ... }
    \item{\%*\%}{\code{signature(x = "diagonalMatrix", y = "dgeMatrix")}: ... }

    \item{solve}{\code{signature(a = "diagonalMatrix", b, ...)}: is
      trivially implemented, of course.}
  }
}
\seealso{
  \code{\linkS4class{ddiMatrix}} and \code{\linkS4class{ldiMatrix}} are
  \dQuote{actual} classes extending \code{"diagonalMatrix"}.
}
\examples{
I5 <- Diagonal(5)
D5 <- Diagonal(x = 10*(1:5))
## trivial (but explicitly defined) methods:
stopifnot(identical(crossprod(I5), I5),
          identical(tcrossprod(I5), I5),
          identical(crossprod(I5, D5), D5),
          identical(tcrossprod(D5, I5), D5),
          identical(solve(D5), solve(D5, I5)),
          all.equal(D5, solve(solve(D5)), tol = 1e-12)
          )
solve(D5)
}
\keyword{classes}
