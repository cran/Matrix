\name{pdMat-class}
\docType{class}
\alias{pdMat-class}
\alias{coef,pdMat-method}
\alias{coerce,pdMat,pdfactor-method}
\alias{coerce,pdMat,matrix-method}
\alias{coerce<-,pdMat,matrix-method}
\alias{coerce<-,pdMat,pdfactor-method}
\alias{dim,pdMat-method}
\alias{formula,pdMat-method}
\alias{names,pdMat-method}
\alias{names<-,pdMat-method}
\alias{show,pdMat-method}
\alias{solve,pdMat-method}
\alias{solve,pdMat,missing-method}
\alias{summary,pdMat-method}
\title{Class pdMat, positive-definite matrices}
\description{A virtual class of parameterized positive-definite
  symmetric matrices.  This class describes the slots and methods that
  actual classes of positive-definite matrices are expected to
  incorporate.  Some classes that inherit from \code{pdMat} have
  additional slots and methods.}
\section{Objects from the Class}{
  Objects of class \code{pdMat} are not constructed directly; only
  objects from classes that inherit from \code{pdMat} are constructed.
}
\section{Slots}{
  \describe{
    \item{\code{form}:}{Object of class \code{"formula"}, holds the
      formula for the object }
    \item{\code{Names}:}{Object of class \code{"character"} holding the
      names of the rows (and columns) of the positive-definite symmetric
      matrix represented by the object.}
    \item{\code{param}:}{\code{"numeric"} - the parameter vector.}
    \item{\code{Ncol}:}{\code{"integer"} - the number of columns (and
      rows) in the matrix.}
    \item{\code{factor}:}{\code{"matrix"} - a square-root factor of the
      matrix.}
    \item{\code{logDet}:}{\code{"numeric"} - the logarithm of the
      determinant of the factor.}
  }
}

\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "pdMat", to = "pdfactor")}:
      extract a square-root factor of the matrix represented by the
      object.  This factor has a \code{logDet} slot giving the
      logarithm of its determinant.  In the case of \code{pdLogChol}
      these are both scalars and the \code{logDet} attribute is the
      logarithm of the absolute value of the factor.}
    \item{corMatrix}{\code{signature(object = "pdMat")}: Extract the
      correlation matrix corresponding to the positive-definite matrix
      represented by the object.  This method is present for back
      compatibility only.  The preferred way of extracting the
      correlation matrix is to coerce the object to the
      \code{"corrmatrix"} class.}
    \item{dim}{\code{signature(x = "pdMat")}: the dimensions of the
      positive-definite matrix represented by the object.} 
    \item{formula}{\code{signature(x = "pdMat")}: extract the formula }
    \item{isInitialized}{\code{signature(object = "pdMat")}:
      \code{TRUE} if the object has been initialized, otherwise \code{FALSE}.}
    \item{logDet}{\code{signature(object = "pdMat", covariate =
	"missing")}: the logarithm of the determinant of the factor of
      the positive-definite matrix represented by the object.} 
    \item{names}{\code{signature(x = "pdMat")}: extract a vector of
      names, which are both the column names and the row names of the
      positive-definite matrix represented by the object.}
    \item{names<-}{\code{signature(x = "pdMat")}: assign the names,}
    \item{pdFactor}{\code{signature(object = "pdMat")}: Extract the
      square root factor of positive-definite symmetric matrix
      represented by the object.  This method is present for back
      compatibility only.  The preferred way of extracting the factor is
      to coerce the object to the \code{"pdfactor"} class.}
    \item{pdMatrix}{\code{signature(object = "pdMat")}: Extract the
      positive-definite symmetric matrix represented by the object.
      This method is present for back compatibility only.  The preferred
      way of extracting the positive-definite symmetric matrix is to
      coerce the object to the \code{"pdmatrix"} class.}
    \item{show}{\code{signature(x = "pdMat")}: show the object.}
    \item{solve}{\code{signature(a = "pdMat", b = "missing")}: Create an
      object of the same class representing the inverse of the
      positive-definite matrix.}
    \item{summary}{\code{signature(object = "pdMat")}:}
  }
}
\references{}
%\author{Douglas Bates \email{bates@stat.wisc.edu}
%  and Saikat DebRoy \email{saikat@stat.wisc.edu}}
\note{
  \code{pdMat} objects are primarily used to represent the variance-covariance
  matrix or the precision matrix of random-effects terms in
  mixed-effects models.  Frequently they are constructed from a formula
  only in the call to the mixed-effects modelling function then assigned
  a value as part of the initialization of the model.
}
\seealso{\code{\link{pdCompSymm-class}}, \code{\link{pdDiag-class}},
  \code{\link{pdLogChol-class}}, \code{\link{pdIdent-class}}
}
%\examples{}
\keyword{classes}
