\name{expand.lu.Matrix.object}
\alias{expand.lu.Matrix.object}
\title{
Expanded LU Decomposition Object
}
\description{
These are objects of class
"expand.lu.Matrix"
They represent the LU decomposition of a matrix.
}
\section{GENERATION}{
This class of objects is returned from the \code{expand.lu.Matrix} function
to represent an LU decomposition.
}
\section{STRUCTURE}{
The following components must be included in a legitimate \code{"expand.lu.Matrix"} object:
\item{l}{
The unit lower triangular/trapezoidal factor.
}
\item{u}{
The upper triangular/trapezoidal factor.
}
\item{permutation}{
The row permutation applied to the original matrix in order to obtain the
factorization.
}}
\seealso{
\code{\link{expand}}, \code{\link{lu.Matrix}}.
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
