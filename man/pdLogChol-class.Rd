\name{pdLogChol-class}
\docType{class}
\alias{pdLogChol-class}
\alias{pdLogChol}
\alias{solve,pdLogChol,missing-method}
\alias{pdLogChol,formula,missing,missing,missing-method}
\title{Class "pdLogChol", positive-definite matrices}
\description{A class of general, positive-definite
  symmetric matrices parameterized by the non-zero elements in the
  Cholesky decomposition.  The diagonal elements are represented by
  their logarithms in the first \code{q} positions of the parameter
  vector.  The strict upper triangle of the factor is in the last
  \code{q(q-1)/2} positions.}
\section{Objects from the Class}{
  Objects of class \code{pdLogChol} can be created by calls of the form
  \code{new("pdLogChol", ...)} or by the generic constructor function
  \code{pdLogChol}.  Frequently the constructor is given a formula only,
  creating an uninitialized \code{pdLogChol} object which is later
  assigned a value.

  \code{pdLogChol} objects are primarily used to represent the
  variance-covariance matrix or the precision matrix of random-effects
  terms in mixed-effects models.
}
\section{Slots}{
  \describe{
    \item{\code{form}:}{Object of class \code{"formula", from class
	"pdMat"}, a formula for the object}
    \item{\code{Names}:}{Object of class \code{"character", from class
	"pdMat"}, names for the rows (and columns) of the
      positive-definite matrix.}
    \item{\code{param}:}{Object of class \code{"numeric", from class
	"pdMat"}, a parameter vector of length
      \eqn{[q(q+1)]/2}{[q(q+1)]/2} where \eqn{q}{q} is \code{Ncol}, the
      number of columns (and rows) in the positive-definite matrix.}
    \item{\code{Ncol}:}{Object of class \code{"integer", from class
	"pdMat"}, number of columns (and rows) in the positive-definite
      matrix.}
    \item{\code{factor}:}{Object of class \code{"matrix", from class
	"pdMat"}, a square root factor of the positive-definite matrix.}
    \item{\code{logDet}:}{Object of class \code{"numeric", from class
	"pdMat"} the logarithm of the absolute value of the determinant
      of the square root factor or, equivalently, half the logarithm of
      the determinant of the positive-definite matrix.}
  }
}
\section{Extends}{
Class \code{"pdMat"}, directly.
}
\section{Methods}{
  \describe{
    \item{EMupdate<-}{\code{signature(x = "pdLogChol", nlev = "numeric",
	value = "matrix")}: update the \code{pdLogChol} object in the EM
      algorithm for a mixed-effects model.}
    \item{LMEgradient}{\code{signature(x = "pdLogChol", A = "matrix",
	nlev = "numeric")}: evaluate the gradient of the log-likelihood
      in a linear mixed-effects model.}
    \item{coef<-}{\code{signature(object = "pdLogChol", value =
	"numeric")}: assign the parameter.}
    \item{coerce}{\code{signature(from = "pdLogChol", to = "pdmatrix")}:
      extract the positive-definite matrix represented by the object.}
    \item{pdgradient}{\code{signature(x = "pdLogChol")}: the gradient of
      the positive definite matrix with respect to the parameter vector.}
    \item{solve}{\code{signature(a = "pdLogChol", b = "missing")}: a
      \code{pdLogChol} object representing the inverse of the
      positive-definite matrix represented by this object.}
    \item{summary}{\code{signature(object = "pdLogChol")}: summarize the
      object.}
  }
}
%\references{}
%\author{Douglas Bates \email{bates@stat.wisc.edu}
%  and Saikat DebRoy \email{saikat@stat.wisc.edu}}
\seealso{\code{\link{pdMat-class}}}
\examples{
m1 <- pdLogChol(~ age)
coef(m1) <- rnorm(3)
print(m1)
solve(m1)
}
\keyword{classes}
