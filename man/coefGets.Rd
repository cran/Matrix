\name{coefGets}
\docType{genericFunction}
\alias{coef<-}
\alias{coef<--methods}
\alias{coef<-,ssclme,numeric-method}
\title{Assign parameters}
\description{
  A generic function to assign the parameter vector in an object.  For
  historical reasons parameters are accessed with \code{coef} and
  assigned with \code{coef<-}.
}
\usage{
coef(object, \dots) <- value
}
\arguments{
  \item{object}{An object that has a numeric parameter vector accessed
    by \code{coef}.}
  \item{\dots}{Additional arguments used by some methods.}
  \item{value}{A numeric vector of the same length as \code{coef(object)}.}
}
%\details{}
\value{
  \code{object} with an updated parameter vector
}
%\references{}
%\author{Douglas Bates \email{bates@stat.wisc.edu}
%  and Saikat DebRoy \email{saikat@stat.wisc.edu}}
%\note{}
%\seealso{}
%\examples{}
\keyword{manip}
