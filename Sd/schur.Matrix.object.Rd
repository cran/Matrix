\name{schur.Matrix.object}
\alias{schur.Matrix.object}
\title{
Schur Decomposition Object
}
\description{
These are objects of class
"schur.Matrix"
They represent the Schur decomposition of a matrix.
}
\section{GENERATION}{
This class of objects is returned from the \code{schur.Matrix} function
to represent an Schur factorization.
}
\section{METHODS}{
The \code{"schur.Matrix"} class of objects has methods for the following generic
functions:
\code{det}, \code{solve}
}
\section{STRUCTURE}{
The following components must be included in a legitimate \code{"schur.Matrix"} 
object:
\item{form}{
a matrix containing either the Schur form of a complex matrix or else the
real Schur form of a non-complex (numeric) matrix. The Schur form of a complex
matrix is a triangular matrix in which the diagonal entries are the 
eigenvalues of the original matrix. The real Schur form of a numeric matrix is
a quasi-triangular matrix in which the diagonal consists of a 1 by 1 block
for each real eigenvalue and a 2 by 2 block for each complex eigenvalue of
the original matrix.
}
\item{vectors}{
the Schur vectors of the matrix, as columns of an orthonormal (or unitary if
complex) matrix.
}}
\seealso{
\code{\link{schur.Matrix}}
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
