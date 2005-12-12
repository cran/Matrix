\name{denseMatrix-class}
\docType{class}
\title{Virtual Class "denseMatrix" of All Dense Matrices}
\alias{denseMatrix-class}
\alias{show,denseMatrix-method}
\alias{coerce,denseMatrix,CsparseMatrix-method}
%
\description{This is the virtual class of all dense (S4) matrices.
  It is the direct superclass of
  \code{\link{ddenseMatrix-class}},
  \code{\link{ldenseMatrix-class}}
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
  see \code{\link{Matrix-class}}.
}
%
\section{Methods}{
  Use \code{\link{showMethods}(class = "denseMatrix", where =
    "package:Matrix")} for an overview of methods.

  Extraction (\code{"["}) methods,
  see \code{\link{[-methods}}.%-> ./Xtrct-methods.Rd
}
\seealso{
  Its superclass \code{\link{Matrix-class}},
  and main subclasses,
  \code{\link{ddenseMatrix-class}}
  \code{\link{sparseMatrix-class}}
}
\examples{
showClass("denseMatrix")
}
\keyword{classes}
