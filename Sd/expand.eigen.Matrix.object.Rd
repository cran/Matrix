\name{expand.eigen.Matrix.object}
\alias{expand.eigen.Matrix.object}
\title{
Expanded Eigenvalue Decomposition Object
}
\description{
These are objects of class
"expand.eigen.Matrix"
They represent the eigenvalue decomposition of a matrix.
}
\section{GENERATION}{
This class of objects is returned from the \code{expand.eigen.Matrix} function
to represent an eigenvalue decomposition.
}
\section{STRUCTURE}{
The following components must be included in a legitimate 
\code{"expand.eigen.Matrix"} object:
\item{values}{
a vector containing the eigenvalues of the matrix.
}
\item{vectors}{
the eigen vectors of the matrix.
A list with two components labeled \code{left} and \code{right}, each of which
contains the corresponding matrix of eigenvectors stored columnwise.
}}
\seealso{
\code{\link{eigen.Matrix}}, \code{\link{expand}}.
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
