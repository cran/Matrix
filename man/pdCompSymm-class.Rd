\name{pdCompSymm-class}
\docType{class}
\alias{pdCompSymm-class}
\alias{pdCompSymm}
\alias{pdCompSymm,formula,missing,missing,missing-method}
\title{Class "pdCompSymm"}
\description{A class of parameterized positive-definite symmetric
  matrices that have compound symmetry.  That is, the diagonal entries
  are constant and the off-diagonal entries are constant.}
\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{new("pdCompSymm", ...)} or
  by calls to the constructor \code{pdCompSymm(...)}
}
\section{Slots}{
  \describe{
    \item{\code{param}:}{Object of class \code{"numeric"} - a parameter
      vector of length 2}
    \item{\code{ncol}:}{Object of class \code{"integer"} - the number of
      columns (and rows) in the matrix.}
    \item{\code{form}:}{Object of class \code{"formula", from class
	"pdMat"} - an optional formula relating the object to a model matrix. }
    \item{\code{Names}:}{Object of class \code{"character", from class
	"pdMat"} - dimnames for the matrix}
  }
}
\section{Extends}{
Class \code{"pdMat"}, directly.
}
\section{Methods}{
  \describe{
    \item{coef}{\code{signature(x = "pdCompSymm")}: get the parameter vector }
    \item{coef<-}{\code{signature(object = "pdCompSymm", value =
	"numeric")}: set the parameter vector}
    \item{coef<-}{\code{signature(object = "pdCompSymm", value =
	"matrix")}: set the parameter vector from a model matrix}
    \item{coerce}{\code{signature(from = "pdCompSymm", to = "matrix")}:
      get the positive-definite symmetric matrix represented by the object.}
    \item{corMatrix}{\code{signature(object = "pdCompSymm")}: get the
      matrix represented by the object as a correlation matrix (with
      standard deviations)}
    \item{dim}{\code{signature(x = "pdCompSymm")}: get the dimensions of
      the matrix represented by the object.}
    \item{isInitialized}{\code{signature(object = "pdCompSymm")}: check
      if the object has been initialized.}
    \item{logDet}{\code{signature(object = "pdCompSymm", covariate =
	"missing")}: get the logarithm of the determinant of the matrix}
    \item{pdFactor}{\code{signature(object = "pdCompSymm")}: get a
      square root factor of the matrix represented by the object.}
    \item{pdMatrix}{\code{signature(object = "pdCompSymm", factor =
	"logical")}: get the positive-definite symmetric matrix
	represented by the object.}
    \item{summary}{\code{signature(object = "pdCompSymm")}:}
  }
}
\author{Douglas Bates \email{bates@stat.wisc.edu}
  and Saikat DebRoy \email{saikat@stat.wisc.edu}}
\seealso{}
\examples{
}
\keyword{classes}
