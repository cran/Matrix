\name{as.Matrix}
\alias{as.Matrix}
\title{
Conversion to Matrix Objects
}
\description{
Converts certain S-PLUS objects to objects of class \code{"Matrix"}. 
}
\usage{
as.Matrix(x)
}
\arguments{
\item{x}{
a matrix, vector, or else an array that is either 2-dimensional or has exactly
two non-unit dimensions.
}}
\value{
an object of class \code{"Matrix"}.
}
\seealso{
\code{\link{Matrix}}
}
\examples{
a <- array( rnorm(4), dim = c(2,1,2))
as.Matrix(a)
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
