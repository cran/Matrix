\name{pdMatrix}
\docType{methods}
\docType{genericFunction}
\alias{pdMatrix}
\alias{pdMatrix-methods}
%\alias{pdMatrix,pdMat-method}
\title{Return the positive-definite matrix}
\description{This generic function extracts the positive-definite matrix
  represented by \code{object}.

  This function is deprecated.  Use coersion to the \code{pdmatrix}
  class instead.
}
\usage{
pdMatrix(object)
}
\arguments{
 \item{object}{an object inheriting from class \code{pdMat}, representing
   a positive definite matrix, which must have been initialized
   (i.e. \code{length(coef(object)) > 0}).} 
}
\value{
  A \code{pdmatrix} object.
}
\seealso{\code{\link{pdmatrix-class}}}
%\examples{}
\keyword{models}
