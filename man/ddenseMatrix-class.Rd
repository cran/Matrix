\name{ddenseMatrix-class}
\docType{class}
\alias{ddenseMatrix-class}
\title{Virtual Class "ddenseMatrix" of Numeric Dense Matrices}
% now separated from ./dgeMatrix-class.Rd - on purpose
% METHODS:
% undocu (FIXME)
\alias{show,ddenseMatrix-method}
\alias{coerce,ddenseMatrix,matrix-method}
%
\alias{cbind2,ddenseMatrix,ddenseMatrix-method}
\alias{cbind2,ddenseMatrix,matrix-method}
\alias{cbind2,matrix,ddenseMatrix-method}
\alias{cbind2,ddenseMatrix,numeric-method}
\alias{cbind2,numeric,ddenseMatrix-method}
\alias{crossprod,ddenseMatrix,missing-method}
\alias{diag,ddenseMatrix-method}
\alias{determinant,ddenseMatrix,missing-method}
\alias{determinant,ddenseMatrix,logical-method}
\alias{lu,ddenseMatrix-method}
\alias{norm,ddenseMatrix,missing-method}
\alias{norm,ddenseMatrix,character-method}
\alias{rcond,ddenseMatrix,missing-method}
\alias{rcond,ddenseMatrix,character-method}
\alias{solve,ddenseMatrix,missing-method}
\alias{solve,ddenseMatrix,ANY-method}
\alias{t,ddenseMatrix-method}
\alias{tcrossprod,ddenseMatrix-method}
%
\description{This is the virtual class of all dense numeric (i.e.,
  \bold{d}ouble, hence \emph{\dQuote{ddense}}) S4 matrices.

  Its most important subclass is the \code{\link{dgeMatrix-class}}.
  %% and now say what the difference is ! __ FIXME __
}
\section{Extends}{
  Class \code{"dMatrix"} directly;
  class \code{"Matrix"}, by the above.
}
\section{Slots}{
  the same slots at its subclass \code{\link{dgeMatrix-class}}, see
  there.
}
\section{Methods}{
  Most methods are implemented via \code{as(*, "dgeMatrix")} and are
  mainly used as \dQuote{fallbacks} when the subclass doesn't need its
  own specialized method.

  Use \code{\link{showMethods}(class = "ddenseMatrix", where =
      "package:Matrix")} for an overview.
}
%\references{}
%\note{}
\seealso{
  \code{\link{Matrix-class}},
  \code{\link{dMatrix-class}},
  \code{\link{dsparseMatrix-class}}
}
\examples{
showClass("ddenseMatrix")

showMethods(class = "ddenseMatrix", where = "package:Matrix")
}
\keyword{classes}
