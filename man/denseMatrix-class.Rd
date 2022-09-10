\name{denseMatrix-class}
\title{Virtual Class "denseMatrix" of All Dense Matrices}
\docType{class}
%
\alias{denseMatrix-class}
%
\alias{-,denseMatrix,missing-method}
\alias{Math,denseMatrix-method}
\alias{as.logical,denseMatrix-method}
\alias{as.numeric,denseMatrix-method}
\alias{as.vector,denseMatrix-method}
\alias{cbind2,denseMatrix,denseMatrix-method}
\alias{cbind2,denseMatrix,matrix-method}
\alias{cbind2,denseMatrix,numeric-method}
%\alias{cbind2,denseMatrix,sparseMatrix-method} % needed in ./cBind.Rd
\alias{cbind2,matrix,denseMatrix-method}
\alias{cbind2,numeric,denseMatrix-method}
\alias{coerce,ANY,denseMatrix-method}
\alias{coerce,denseMatrix,CsparseMatrix-method}
\alias{coerce,denseMatrix,RsparseMatrix-method}
\alias{coerce,denseMatrix,TsparseMatrix-method}
\alias{coerce,denseMatrix,dMatrix-method}
\alias{coerce,denseMatrix,ddenseMatrix-method}
\alias{coerce,denseMatrix,dsparseMatrix-method}
\alias{coerce,denseMatrix,generalMatrix-method}
\alias{coerce,denseMatrix,lMatrix-method}
\alias{coerce,denseMatrix,ldenseMatrix-method}
\alias{coerce,denseMatrix,lsparseMatrix-method}
\alias{coerce,denseMatrix,matrix-method}
\alias{coerce,denseMatrix,nMatrix-method}
\alias{coerce,denseMatrix,ndenseMatrix-method}
\alias{coerce,denseMatrix,nsparseMatrix-method}
\alias{coerce,denseMatrix,packedMatrix-method}
\alias{coerce,denseMatrix,sparseMatrix-method}
\alias{coerce,denseMatrix,unpackedMatrix-method}
\alias{coerce,denseMatrix,vector-method}
\alias{coerce,matrix,denseMatrix-method}
\alias{coerce,numLike,denseMatrix-method}
\alias{dim<-,denseMatrix-method}
\alias{log,denseMatrix-method}
\alias{rbind2,denseMatrix,denseMatrix-method}
\alias{rbind2,denseMatrix,matrix-method}
\alias{rbind2,denseMatrix,numeric-method}
%\alias{rbind2,denseMatrix,sparseMatrix-method} % needed in ./cBind.Rd
\alias{rbind2,matrix,denseMatrix-method}
\alias{rbind2,numeric,denseMatrix-method}
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
\keyword{classes}
