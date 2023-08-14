\name{denseMatrix-class}
\title{Virtual Class "denseMatrix" of All Dense Matrices}
%
\docType{class}
\keyword{array}
\keyword{classes}
%
\alias{denseMatrix-class}
%
\alias{-,denseMatrix,missing-method}
\alias{Math,denseMatrix-method}
\alias{coerce,ANY,denseMatrix-method}
\alias{coerce,matrix,denseMatrix-method}
\alias{coerce,vector,denseMatrix-method}
\alias{dim<-,denseMatrix-method}
\alias{log,denseMatrix-method}
\alias{mean,denseMatrix-method}
\alias{rep,denseMatrix-method}
\alias{show,denseMatrix-method}
%
\description{This is the virtual class of all dense (S4) matrices.
  It partitions into two subclasses
  \code{\linkS4class{packedMatrix}} and
  \code{\linkS4class{unpackedMatrix}}.
  Alternatively into the (currently) three subclasses
  \code{\linkS4class{ddenseMatrix}},
  \code{\linkS4class{ldenseMatrix}}, and
  \code{\linkS4class{ndenseMatrix}}.

  \code{denseMatrix} is (hence) the direct superclass of these (\eqn{2+3 = 5}) classes.
}
\section{Extends}{
  class \code{"Matrix"} directly.
}
\section{Slots}{
  exactly those of its superclass \code{"\linkS4class{Matrix}"}, i.e.,
  \code{"Dim"} and \code{"Dimnames"}.
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
