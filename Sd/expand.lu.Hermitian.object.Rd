\name{expand.lu.Hermitian.object}
\alias{expand.lu.Hermitian.object}
\title{
Expanded Symmetric-Indefinite Decomposition Object
}
\description{
These are objects of class
"expand.lu.Hermitian"
They represent the symmetric indefinite decomposition 
of a real symmetric or complex Hermitian matrix.
}
\section{GENERATION}{
This class of objects is returned from the 
\code{expand.lu.Hermitian} function to represent a symmetric 
indefinite decomposition.
}
\section{STRUCTURE}{
The following components must be included in a legitimate 
\code{"expand.lu.Hermitian"} object:
\item{l}{
The unit lower triangular/trapezoidal factor.
}
\item{u}{
The upper triangular/trapezoidal factor.
}
\item{permutation}{
The row permutation applied symmetrically to the original matrix in order to 
obtain the factorization.
}}
\seealso{
\code{\link{lu.Hermitian}}, \code{\link{expand}}.
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
