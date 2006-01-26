\name{dgCMatrix-class}
\docType{class}
\title{Compressed, sparse, column-oriented numeric matrices}
\alias{dgCMatrix-class}
\alias{\%*\%,dgCMatrix,dgeMatrix-method}
\alias{\%*\%,dgCMatrix,matrix-method}
\alias{\%*\%,dgeMatrix,dgCMatrix-method}
\alias{\%*\%,matrix,dgCMatrix-method}
\alias{coerce,matrix,dgCMatrix-method}
\alias{coerce,dgeMatrix,dgCMatrix-method}
\alias{coerce,dgCMatrix,dgTMatrix-method}
\alias{coerce,dgCMatrix,dgBCMatrix-method}
\alias{coerce,dgCMatrix,lgCMatrix-method}
\alias{coerce,dgCMatrix,matrix-method}
\alias{coerce,dgCMatrix,dgeMatrix-method}
\alias{crossprod,dgCMatrix,missing-method}
\alias{crossprod,dgCMatrix,matrix-method}
\alias{crossprod,dgCMatrix,dgeMatrix-method}
\alias{diag,dgCMatrix-method}
\alias{dim,dgCMatrix-method}
\alias{image,dgCMatrix-method}
\alias{t,dgCMatrix-method}
%% Group methods --------- FIXME: These are not tested yet (or documented)
\alias{Arith,numeric,dgCMatrix-method}
\alias{Arith,dgCMatrix,numeric-method}
%\alias{Math2,dgCMatrix,numeric-method}
\alias{Math,dgCMatrix-method}
% for silly reasons, need these 2+3 as well:
\alias{round,dgCMatrix,numeric-method}
\alias{signif,dgCMatrix,numeric-method}
\alias{log,dgCMatrix-method}
\alias{gamma,dgCMatrix-method}
\alias{lgamma,dgCMatrix-method}
%
\description{The \code{dgCMatrix} class is a class of sparse numeric
  matrices in the compressed, sparse, column-oriented format.  In this
  implementation the non-zero elements in the columns are sorted into
  increasing row order.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("dgCMatrix", ...)} or
  often more easily via \code{\link{Matrix}(*, sparse = TRUE)}.
}
\section{Slots}{
  \describe{
    \item{\code{x}:}{Object of class \code{"numeric"} - the non-zero
      elements of the matrix.}
    \item{\dots}{all other slots are inherited from the superclass
      \code{"CsparseMatrix"}, see \code{\link{CsparseMatrix-class}}.
    }
  }
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "matrix", to = "dgCMatrix")}}
    \item{coerce}{\code{signature(from = "dgCMatrix", to = "matrix")}}
    \item{coerce}{\code{signature(from = "dgCMatrix", to = "dgBCMatrix")}}
    \item{coerce}{\code{signature(from = "dgCMatrix", to = "dgTMatrix")}}
    \item{crossprod}{\code{signature(x = "dgCMatrix", y = "missing")}
      returns \code{t(x) \%*\% x} as an \code{dsCMatrix} object.}
    \item{crossprod}{\code{signature(x = "dgCMatrix", y = "matrix")}: ... }
    \item{crossprod}{\code{signature(x = "dgCMatrix", y = "numeric")}: ... }
    \item{diag}{\code{signature(x = "dgCMatrix")}: returns the diagonal
      of \code{x}}
    \item{dim}{\code{signature(x = "dgCMatrix")}: returns the dimensions
      of \code{x}}
    \item{image}{\code{signature(x = "dgCMatrix")}: plots an image of
      \code{x} using the \code{\link[lattice]{levelplot}} function}
  }
}
%\references{}
%\author{}
%\note{}
\seealso{
  \code{\link{dsCMatrix-class}}, \code{\link{dtCMatrix-class}}
}
\examples{
(m <- Matrix(c(0,0,2:0), 3,5))
str(m)
m[,1]
}
\keyword{classes}
\keyword{algebra}
