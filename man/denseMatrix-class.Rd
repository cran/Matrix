\name{denseMatrix-class}
\docType{class}
\title{Virtual Class "denseMatrix" of All Dense Matrices}
\alias{denseMatrix-class}
%% Group methods
\alias{Math,denseMatrix-method}
%
\alias{show,denseMatrix-method}
\alias{coerce,ANY,denseMatrix-method}
\alias{coerce,denseMatrix,CsparseMatrix-method}
\alias{coerce,denseMatrix,TsparseMatrix-method}
\alias{coerce,denseMatrix,generalMatrix-method}
\alias{coerce,denseMatrix,sparseMatrix-method}
\alias{cbind2,denseMatrix,denseMatrix-method}
\alias{cbind2,denseMatrix,matrix-method}
\alias{cbind2,matrix,denseMatrix-method}
\alias{cbind2,denseMatrix,numeric-method}
\alias{cbind2,numeric,denseMatrix-method}
\alias{rbind2,denseMatrix,denseMatrix-method}
\alias{rbind2,denseMatrix,matrix-method}
\alias{rbind2,matrix,denseMatrix-method}
\alias{rbind2,denseMatrix,numeric-method}
\alias{rbind2,numeric,denseMatrix-method}
\alias{dim<-,denseMatrix-method}
\alias{rcond,denseMatrix,character-method}
%
\description{This is the virtual class of all dense (S4) matrices.
  It is the direct superclass of
  \code{\linkS4class{ddenseMatrix}},
  \code{\linkS4class{ldenseMatrix}}
}
\section{Extends}{
  class \code{"Matrix"} directly.
}
\section{Slots}{
  %% MM :adding the slots here because "R CMD check" gave warning
  %%  _WHY_ was this not necessary in ./ddenseMatrix-class.Rd ?
  %%   or ./sparseMatrix-class.Rd  ???
  exactly those of its superclass \code{"Matrix"}, i.e.,
  \describe{
    \item{\code{factors}:}{Object of class \code{"list"} - a list
      of factorizations of the matrix.}
    \item{\code{Dim}:}{length 2 \code{"integer"}}
    \item{\code{Dimnames}:}{\code{"list"} of length 2,}
  }
  %% not sufficient: \code{Dim}, and \code{Dimnames},
  see \code{\linkS4class{Matrix}}.
}
%
\section{Methods}{
  Use \code{\link{showMethods}(class = "denseMatrix", where =
    "package:Matrix")} for an overview of methods.

  Extraction (\code{"["}) methods,
  see \code{\link{[-methods}}.%-> ./Xtrct-methods.Rd
}
\seealso{
  \code{\link{colSums}}, \code{\link{kronecker}}, and other such methods
  with own help pages.

  Its superclass \code{\linkS4class{Matrix}}, and main subclasses,
  \code{\linkS4class{ddenseMatrix}} and \code{\linkS4class{sparseMatrix}}.
}
\examples{
showClass("denseMatrix")
}
\keyword{classes}
