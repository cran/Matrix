\name{eigen.Matrix.object}
\alias{eigen.Matrix.object}
\title{
Eigenvalue Decomposition Object
}
\description{
These are objects of class
"eigen.Matrix"
They represent the eigenvalue decomposition of a matrix.
}
\section{GENERATION}{
This class of objects is returned from the \code{eigen.Matrix} function
to represent an eigenvalue decomposition.
}
\section{METHODS}{
The \code{"eigen.Matrix"} class of objects has methods for the following generic
functions:
\code{det}, \code{expand}, \code{facmul}, \code{solve}
}
\section{STRUCTURE}{
The following components must be included in a legitimate
\code{"eigen.Matrix"} object:
\item{values}{
a vector containing the eigenvalues of the matrix.
}
\item{vectors}{
the eigenvectors of the matrix, represented as a list with two components 
labeled \code{left} and \code{right}.
If the underlying matrix is complex or real with real eigenvalues,
these components are just the matrices of left and right eigenvectors 
of the underlying matrix stored columnwise.
If the matrix is real but has complex eigenvalues, then these eigenvalues 
will appear in complex conjugate pairs in \code{values}. Rather than store
a complex matrix of eigenvectors, the real and imaginary parts of the 
eigenvectors for complex conjugate pairs are stored in the corresponding 
columns of \code{left} and \code{right}.
}}
\seealso{
\code{\link{eigen.Hermitian}}, \code{\link{eigen.Matrix}}.
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
