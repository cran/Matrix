\name{lu.Hermitian.object}
\alias{lu.Hermitian.object}
\title{
Symmetric Indefinite Factorization Object
}
\description{
These are objects of class
"lu.Hermitian"
They represent the symmetric indefinite factorization of a 
real symmetric or complex Hermitian matrix.
}
\section{GENERATION}{
This class of objects is returned from the \code{lu.Hermitian}
function to represent a symmetric indefinite factorization.
}
\section{METHODS}{
The \code{"lu.Hermitian"} class of objects has methods for the 
following generic functions:
\code{det}, \code{expand}, \code{facmul}, \code{rcond}, \code{solve}
}
\section{STRUCTURE}{
The following component must be included in a legitimate \code{"lu.Hermitian"} 
object:
\item{factors}{
a matrix in which the triangular and block diagonal factors are stored.
There is an attribute \code{"uplo"} associated with this matrix indicating
whether or not the triangular factor is upper or lower triangular.
}
\item{pivot}{
the permutation applied to the original matrix to achieve the
symmetric indefinite factorization.
}}
\seealso{
\code{\link{lu.Hermitian}}
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
