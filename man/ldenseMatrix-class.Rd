\name{ldenseMatrix-class}
\title{Virtual Class "ldenseMatrix" of Dense Logical Matrices}
%
\docType{class}
\keyword{array}
\keyword{classes}
%
\alias{ldenseMatrix-class}
%
\alias{!,ldenseMatrix-method}
\alias{&,ldenseMatrix,ddiMatrix-method}
\alias{&,ldenseMatrix,ldiMatrix-method}
\alias{&,ldenseMatrix,ndiMatrix-method}
\alias{*,ldenseMatrix,ddiMatrix-method}
\alias{*,ldenseMatrix,ldiMatrix-method}
\alias{*,ldenseMatrix,ndiMatrix-method}
\alias{-,ldenseMatrix,missing-method}
\alias{-,ndenseMatrix,missing-method}
\alias{Logic,ldenseMatrix,lsparseMatrix-method}
\alias{Ops,ldenseMatrix,ldenseMatrix-method}
\alias{^,ldenseMatrix,ddiMatrix-method}
\alias{^,ldenseMatrix,ldiMatrix-method}
\alias{^,ldenseMatrix,ndiMatrix-method}
\alias{coerce,matrix,ldenseMatrix-method}
\alias{coerce,vector,ldenseMatrix-method}
\alias{which,ldenseMatrix-method}
%
\description{
  \code{ldenseMatrix} is the virtual class of all dense \bold{l}ogical
  (S4) matrices.  It extends both \code{\linkS4class{denseMatrix}}
  and \code{\linkS4class{lMatrix}} directly.
}
\section{Slots}{
  \describe{
    \item{\code{x}:}{logical vector containing the entries of the matrix.}
    \item{\code{Dim}, \code{Dimnames}:}{see \code{\linkS4class{Matrix}}.}
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
    \item{as.vector}{\code{signature(x = "ldenseMatrix", mode = "missing")}: ...}
    \item{which}{\code{signature(x = "ndenseMatrix")}, semantically
      equivalent to \pkg{base} function \code{\link{which}(x, arr.ind)};
      for details, see the \code{\linkS4class{lMatrix}} class documentation.}
  }
}
\seealso{
 Class \code{\linkS4class{lgeMatrix}} and the other subclasses.
}
\examples{
showClass("ldenseMatrix")

as(diag(3) > 0, "ldenseMatrix")
}
