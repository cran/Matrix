\name{eigen.Hermitian.object}
\alias{eigen.Hermitian.object}
\title{
Hermitian Eigenvalue Decomposition Object
}
\description{
These are objects of class
"eigen.Hermitian"
They represent the eigenvalue decomposition of a symmetric or Hermitian matrix.
}
\section{GENERATION}{
This class of objects is returned from the 
\code{eigen.Hermitian} function to represent a Hermitian eigenvalue decomposition.
}
\section{METHODS}{
The \code{"eigen.Hermitian"} classes of objects have 
methods for the following generic functions:
\code{det}, \code{norm}, \code{rcond}, \code{solve}
}
\section{STRUCTURE}{
The following components must be included in a legitimate \code{"eigen.Hermitian"} 
object:
\item{values}{
a vector containing the eigenvalues of the matrix.
}
\item{vectors}{
the orthonormal matrix of eigenvectors, stored columnwise.
}}
\seealso{
\code{\link{eigen.Matrix}}.
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
