\name{lu.Matrix.object}
\alias{lu.Matrix.object}
\title{
LU Decomposition Object
}
\description{
These are objects of class
"lu.Matrix"
They represent the LU decomposition of a Matrix.
}
\section{GENERATION}{
This class of objects is returned from the \code{lu.Matrix} function
to represent an LU decomposition.
}
\section{METHODS}{
The \code{"lu.Matrix"} class of objects has methods for the following generic
functions:
\code{det}, \code{expand}, \code{facmul}, \code{rcond}, \code{solve}
}
\section{STRUCTURE}{
The following component must be included in a legitimate \code{"lu.Matrix"} object:
\item{factors}{
a matrix in which the triangular factors are stored.
}
\item{pivot}{
the row permutation applied to the original matrix to achieve LU decomposition.
}}
\seealso{
\code{\link{lu.Matrix}}
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
