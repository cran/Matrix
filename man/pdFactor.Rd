\name{pdFactor}
\docType{methods}
\docType{genericFunction}
\alias{pdFactor}
\alias{pdFactor-methods}
%\alias{pdFactor,pdMat-method}
\title{Square-root factor}
\description{This generic function extracts the square-root factor of
  the positive-definite matrix represented by \code{object}. Letting
  \eqn{\Sigma}{S} denote a positive-definite matrix, a square-root
  factor of \eqn{\Sigma}{S} is any square matrix \eqn{L}{L} such that
  \eqn{\Sigma = L'L}{S = L'L}. This function extracts \eqn{L}.

  This function is deprecated.  Use coersion to the \code{pdfactor}
  class instead.
}
\usage{
pdFactor(object)
}
\arguments{
 \item{object}{an object inheriting from class \code{pdMat}, representing
   a positive definite matrix, which must have been initialized
   (i.e. \code{length(coef(object)) > 0}).} 
}
\value{
  A \code{pdfactor} object.
}
\seealso{\code{\link{pdfactor-class}}}
%\examples{}
\keyword{models}
