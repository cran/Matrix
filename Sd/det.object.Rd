\name{det.object}
\alias{det.object}
\title{
Determinant Object
}
\description{
These are objects of class
"det"
They represent the determinant of a Matrix.
}
\section{GENERATION}{
This class of objects is returned from \code{det} methods.
}
\section{METHODS}{
The \code{"det"} class of objects has methods for the following generic functions:
\code{print}
}
\section{STRUCTURE}{
The following components must be included in a legitimate \code{"det"} object:
\item{modulus}{
Either the modulus of the determinant of \code{x} or its logarithm. 
A logical attribute \code{"logarithm"} is included to indicate whether the true
modulus or its logarithm is given.
}
\item{sign}{
the sign of the determinant of \code{x}. For complex matrices, this will 
usually be a complex number with unit modulus.


A copy of the call to \code{det} is returned as an attribute, as well as 
workspace information if provided by the method.
}}
\seealso{
\code{\link{det}}, \code{\link{sign}}.
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
