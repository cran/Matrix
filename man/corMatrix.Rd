\name{corMatrix}
\docType{genericFunction}
\alias{corMatrix}
\alias{corMatrix,pdMat-method}
\title{Generate the correlation matrix}
\description{
  This generic function generates a correlation matrix from an object of
  a suitable class.

  This function is deprecated.  Use coersion to the \code{corrmatrix}
  class instead.
}
\usage{
corMatrix(object, ...)
} 
\arguments{
  \item{object}{An object of a suitable class}
  \item{\dots}{Some methods for this generic take additional, optional
    arguments.}
}
%\details{}
\value{
  An object of class \code{\link[lme4:corrmatrix-class]{corrmatrix}}.
}
%\references{}
%\author{Douglas Bates \email{bates@stat.wisc.edu} and Saikat DebRoy \email{saikat@stat.wisc.edu}}
%\note{}
\seealso{\code{\link{corrmatrix-class}}}
%\examples{}
\keyword{models}
