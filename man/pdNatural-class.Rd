\name{pdNatural-class}
\docType{class}
\alias{pdNatural-class}
\alias{pdNatural}
\alias{pdNatural,formula,missing,missing,missing-method}
\alias{pdNatural,pdMat,missing,missing,missing-method}
\alias{solve,pdNatural,missing-method}
\title{Class "pdNatural", positive-definite matrices}
\description{A class of general, positive-definite
  symmetric matrices parameterized by the logarithm of the diagonal
  elements and Fisher's z transformation of the correlations.}
\section{Objects from the Class}{

  Objects can be created by calls of the form \code{new("pdNatural",
    ...)} or by the generic constructor function \code{pdNatural}.
  Frequently the constructor is given a formula only, creating an
  uninitialized \code{pdNatural} object which is later assigned a
  value.

  \code{pdNatural} objects are primarily used to represent the
  variance-covariance matrix or the precision matrix of random-effects
  terms in mixed-effects models.
}
\section{Slots}{
  \describe{
    \item{\code{param}:}{Object of class \code{"numeric", from class
	"pdMat"}, a parameter vector of length
      \eqn{[q(q+1)]/2}{[q(q+1)]/2} where \eqn{q}{q} is the number of
      rows (and columns) in the positive-definite matrix.}
    \item{\code{form}:}{Object of class \code{"formula", from class
	"pdMat"}, a formula for the object}
    \item{\code{Names}:}{Object of class \code{"character", from class
	"pdMat"}, names for the rows (and columns) of the
      positive-definite matrix.}
  }
}
\section{Extends}{
Class \code{"pdMat"}, directly.
}
\section{Methods}{
  \describe{
    \item{EMupdate<-}{\code{signature(x = "pdNatural", nlev = "numeric",
	value = "matrix")}: update the \code{pdNatural} object in the EM
      algorithm for a mixed-effects model.}
    \item{LMEgradient}{\code{signature(x = "pdNatural", A = "matrix",
	nlev = "numeric")}: evaluate the gradient of the log-likelihood
      in a linear mixed-effects model.}
    \item{coef}{\code{signature(object = "pdNatural")}: extract the
      parameter.}
    \item{coef<-}{\code{signature(object = "pdNatural", value =
	"numeric")}: assign the parameter.}
    \item{coerce}{\code{signature(from = "pdNatural", to = "pdmatrix")}:
      extract the positive-definite matrix represented by the object.}
    \item{coerce}{\code{signature(from = "pdNatural", to = "pdfactor")}:
      extract a square-root factor of the matrix represented by the
      object.  This factor has a \code{logDet} attribute giving the
      logarithm of its determinant.  In the case of \code{pdNatural}
      these are both scalars and the \code{logDet} attribute is the
      logarithm of the absolute value of the factor.}
    \item{dim}{\code{signature(x = "pdNatural")}: the dimensions of the
      positive-definite matrix represented by the object.} 
    \item{isInitialized}{\code{signature(object = "pdNatural")}: a
      logical scalar indicating if the object is initialized. }
    \item{logDet}{\code{signature(object = "pdNatural", covariate =
	"missing")}: the logarithm of the determinant of the factor of
      the positive-definite matrix represented by the object.} 
    \item{pdgradient}{\code{signature(x = "pdNatural")}: the gradient of
      the positive definite matrix with respect to the parameter vector.}
    \item{solve}{\code{signature(a = "pdNatural", b = "missing")}: a
      \code{pdNatural} object representing the inverse of the
      positive-definite matrix represented by this object.}
    \item{summary}{\code{signature(object = "pdNatural")}: summarize the
      object.}
  }
}
%\references{}
%\author{Douglas Bates \email{bates@stat.wisc.edu}
%  and Saikat DebRoy \email{saikat@stat.wisc.edu}}
\seealso{\code{\link{pdMat-class}}}
%\examples{}
\keyword{classes}
