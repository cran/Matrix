\name{pdIdent-class}
\docType{class}
\alias{pdIdent-class}
\alias{pdIdent}
\alias{solve,pdIdent,missing-method}
\alias{pdIdent,formula,missing,missing,missing-method}
\title{Class "pdIdent"}
\description{A class of parameterized positive-definite
  matrices that are multiples of the identity matrix.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("pdIdent", ...)} or
  by calls to the constructor function \code{pdIdent}.
}
\section{Slots}{
  \describe{
    \item{\code{form}:}{Object of class \code{"formula", from class
	"pdMat"} an optional formula used to generate the object.}
    \item{\code{Names}:}{Object of class \code{"character", from class
	"pdMat"} row (and column) names for the matrix.}
    \item{\code{param}:}{Object of class \code{"numeric", from class
	"pdMat"} the unconstrained, numeric parameter vector of length
      \code{Ncol}.  The elements of the parameter vector are the
      logarithms of the square roots of the diagonal elements of the matrix. }
    \item{\code{Ncol}:}{Object of class \code{"integer", from class
	"pdMat"} the number of rows and columns in the matrix. }
    \item{\code{factor}:}{Object of class \code{"matrix", from class
	"pdMat"} a square-root factor of the matrix - also a diagonal
      matrix.}
    \item{\code{logDet}:}{Object of class \code{"numeric", from class
	"pdMat"} the logarithm of the absolute value of the determinant
      of the factor.}
  }
}
\section{Extends}{
Class \code{"pdMat"}, directly.
}
\section{Methods}{
  \describe{
    \item{coef}{\code{signature(x = "pdIdent")}: extract the parameters.}
    \item{coef<-}{\code{signature(object = "pdIdent", value =
	"numeric")}: assign the parameters.}
    \item{coerce}{\code{signature(from = "pdIdent", to = "matrix")}:
      extract the matrix.}
    \item{corMatrix}{\code{signature(object = "pdIdent")}: extract the
      correlation matrix. }
    \item{dim}{\code{signature(x = "pdIdent")}: extract the dimensions.}
    \item{isInitialized}{\code{signature(object = "pdIdent")}: check if initialized.}
    \item{logDet}{\code{signature(object = "pdIdent", covariate =
	"missing")}: return the logarithm of the determinant of the matrix.}
    \item{solve}{\code{signature(a = "pdIdent", b = "missing")}: invert
      the matrix returning a \code{pdIdent} object.}
  }
}
%\references{}
%\author{Douglas Bates \email{bates@stat.wisc.edu} and Saikat DebRoy \email{saikat@stat.wisc.edu}}
%\note{}
\seealso{\code{\link{pdMat-class}}}
\examples{
}
\keyword{classes}
