\name{ldenseMatrix-class}
\docType{class}
\alias{ldenseMatrix-class}
\alias{!,ldenseMatrix-method}
\alias{coerce,matrix,ldenseMatrix-method}
\alias{coerce,ldenseMatrix,matrix-method}
%%-- These are already in ./lsparseMatrix-classes.Rd :
%\alias{\%*\%,lsparseMatrix,ldenseMatrix-method}
%\alias{\%*\%,ldenseMatrix,lsparseMatrix-method}
%\alias{crossprod,lsparseMatrix,ldenseMatrix-method}
%\alias{crossprod,ldenseMatrix,lsparseMatrix-method}
%
\title{Virtual Class "ldenseMatrix" of Dense Logical Matrices}
\description{
  \code{ldenseMatrix} is the virtual class of all dense \bold{l}ogical
  (S4) matrices.  It extends both \code{\link{denseMatrix-class}}
  and \code{\link{lMatrix-class}} directly.
}
\section{Slots}{
  \describe{
    \item{\code{x}:}{logical vector containing the entries of the matrix.}
    \item{\code{factors}:}{Object of class \code{"list"}, possibly
      containing matrix factorizations.}
    \item{\code{Dim}, \code{Dimnames}:}{see \code{\link{Matrix-class}}.}
  }
}
\section{Extends}{
Class \code{"lMatrix"}, directly.
Class \code{"denseMatrix"}, directly.
Class \code{"Matrix"}, by class \code{"lMatrix"}.
Class \code{"Matrix"}, by class \code{"denseMatrix"}.
}
\section{Methods}{
  \describe{
    \item{\%*\%}{\code{signature(x = "lsparseMatrix", y = "ldenseMatrix")}: ... }
    \item{\%*\%}{\code{signature(x = "ldenseMatrix", y = "lsparseMatrix")}: ... }
    \item{coerce}{\code{signature(from = "matrix", to = "ldenseMatrix")}: ... }
    \item{coerce}{\code{signature(from = "ldenseMatrix", to = "matrix")}: ... }
    \item{crossprod}{\code{signature(x = "lsparseMatrix", y = "ldenseMatrix")}: ... }
    \item{crossprod}{\code{signature(x = "ldenseMatrix", y = "lsparseMatrix")}: ... }
  }
}
\seealso{
 \code{\link{lgeMatrix-class}} and the other subclasses.
}
\examples{
showClass("ldenseMatrix")
}
\keyword{classes}
