\name{svd.Matrix.object}
\alias{svd.Matrix.object}
\title{
Singular Value Decomposition Object
}
\description{
These are objects of class
"svd.Matrix"
They represent the singular value decomposition of a matrix.
}
\section{GENERATION}{
This class of objects is returned from the \code{svd.Matrix} function
to represent a singular value decomposition.
}
\section{METHODS}{
The \code{"svd.Matrix"} class of objects has methods for the following generic
functions:
\code{det}, \code{norm}, \code{rcond}, \code{solve}
}
\section{STRUCTURE}{
The following components must be included in a legitimate \code{"svd.Matrix"} 
object:
\item{values}{
the vector of singular values (in decreasing order of magnitude).
}
\item{vectors}{
the singular vectors of the matrix. 
A list with two components labeled \code{left} and \code{right}, each of which
contains the corresponding orthonormal matrix of singular vectors stored 
columnwise.
}}
\seealso{
\code{\link{svd.Matrix}}.
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
