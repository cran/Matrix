\name{unpack}
\alias{unpack}
\title{
Full Storage Representation for Packed Matrices --- Generic Function
}
\description{
Expands decompositions stored in compact form into matrix factors.
Methods are available for the following classes from the 
\code{Matrix} library: \code{"unpack.Diagonal"}, \code{"unpack.Identity"}, 
\code{"unpack.ColPermutation"}, \code{"unpack.RowPermutation"}
}
\usage{
unpack(x, \dots)
}
\arguments{
\item{x}{
a matrix stored in packed form. 
}
\item{...}{
additional arguments needed (not currently used in any method).
}}
\value{
the full-storage representation of \code{x}.
}
\examples{
library(Matrix)
x <- Diagonal( 1:3)
x
unpack(x)
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
