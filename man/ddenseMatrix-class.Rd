\name{ddenseMatrix-class}
\title{Virtual Class "ddenseMatrix" of Numeric Dense Matrices}
%
\docType{class}
\keyword{array}
\keyword{classes}
%
\alias{ddenseMatrix-class}
%
\alias{&,ddenseMatrix,ddiMatrix-method}
\alias{&,ddenseMatrix,ldiMatrix-method}
\alias{&,ddenseMatrix,ndiMatrix-method}
\alias{*,ddenseMatrix,ddiMatrix-method}
\alias{*,ddenseMatrix,ldiMatrix-method}
\alias{*,ddenseMatrix,ndiMatrix-method}
\alias{Arith,ddenseMatrix,logical-method}
\alias{Arith,ddenseMatrix,numeric-method}
\alias{Arith,ddenseMatrix,sparseVector-method}
\alias{Arith,logical,ddenseMatrix-method}
\alias{Arith,numeric,ddenseMatrix-method}
\alias{^,ddenseMatrix,ddiMatrix-method}
\alias{^,ddenseMatrix,ldiMatrix-method}
\alias{^,ddenseMatrix,ndiMatrix-method}
\alias{coerce,matrix,ddenseMatrix-method}
\alias{coerce,vector,ddenseMatrix-method}
%
\description{This is the virtual class of all dense numeric (i.e.,
  \bold{d}ouble, hence \emph{\dQuote{ddense}}) S4 matrices.

  Its most important subclass is the \code{\linkS4class{dgeMatrix}} class.
  %% and now say what the difference is ! __ FIXME __
}
\section{Extends}{
  Class \code{"dMatrix"} directly;
  class \code{"Matrix"}, by the above.
}
\section{Slots}{
  the same slots at its subclass \code{\linkS4class{dgeMatrix}}, see
  there.
}
\section{Methods}{
  Most methods are implemented via \code{as(*, "generalMatrix")} and are
  mainly used as \dQuote{fallbacks} when the subclass doesn't need its
  own specialized method.

  Use \code{\link{showMethods}(class = "ddenseMatrix", where =
      "package:Matrix")} for an overview.
}
%\references{}
%\note{}
\seealso{
  The virtual classes \code{\linkS4class{Matrix}},
  \code{\linkS4class{dMatrix}}, and \code{\linkS4class{dsparseMatrix}}.
}
\examples{
showClass("ddenseMatrix")

showMethods(class = "ddenseMatrix", where = "package:Matrix")
}
