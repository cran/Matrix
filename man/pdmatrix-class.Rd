\name{pdmatrix-class}
\docType{class}
\alias{pdmatrix-class}
\alias{coerce,matrix,pdmatrix-method}
\alias{coerce,pdCompSymm,pdmatrix-method}
\alias{coerce,pdDiag,pdmatrix-method}
\alias{coerce,pdIdent,pdmatrix-method}
\alias{coerce,pdLogChol,pdmatrix-method}
\alias{coerce,pdMat,pdmatrix-method}
\alias{coerce,pdNatural,pdmatrix-method}
\alias{coerce,pdfactor,pdmatrix-method}
\alias{coerce,pdmatrix,pdfactor-method}
\alias{coerce<-,pdDiag,pdmatrix-method}
\alias{coerce<-,pdIdent,pdmatrix-method}
\alias{coerce<-,pdLogChol,pdmatrix-method}
\alias{coerce<-,pdMat,pdmatrix-method}
\alias{coerce<-,pdNatural,pdmatrix-method}
\alias{coerce<-,pdmatrix,pdfactor-method}
\title{Positive-definite matrices }
\description{The class of positive-definite, symmetric matrices.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("pdmatrix", ...)}.
  More frequently they are created by coercing a \code{pdMat}
  object to the \code{pdmatrix} class.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"matrix"}, which should
      be a numeric matrix.}
  }
}
\section{Extends}{
  Class \code{"matrix"}, from data part.
  Class \code{"structure"}, by class "matrix".
  Class \code{"array"}, by class "matrix".
  Class \code{"vector"}, by class "matrix".
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "pdDiag", to = "pdmatrix")}: extract the positive-definite matrix }
    \item{coerce}{\code{signature(from = "pdIdent", to = "pdmatrix")}:
      extract the positive-definite matrix}
    \item{coerce}{\code{signature(from = "pdLogChol", to = "pdmatrix")}:
      extract the positive-definite matrix}
    \item{coerce}{\code{signature(from = "pdNatural", to = "pdmatrix")}:
      extract the positive-definite matrix}
%    \item{coerce}{\code{signature(from = "pdScalar", to = "pdmatrix")}:
%      extract the positive-definite matrix}
    \item{coerce}{\code{signature(from = "matrix", to = "pdmatrix")}:
      create a \code{pdmatrix} object from a matrix (does not check for
      symmetry or positive-definiteness)}
    \item{coerce}{\code{signature(from = "pdfactor", to = "pdmatrix")}:
      create a positive-definite matrix as the crossproduct of the factor.}
    \item{coerce}{\code{signature(from = "pdmatrix", to = "pdfactor")}:
      convert a \code{pdmatrix} object to its Cholesky factor}
  }
}
%\references{}
%\author{Douglas Bates \email{bates@stat.wisc.edu} and Saikat DebRoy \email{saikat@stat.wisc.edu}}
%\note{}
\seealso{\code{\link{pdfactor-class}}, \code{\link{pdMat-class}}}
%\examples{}
\keyword{classes}
