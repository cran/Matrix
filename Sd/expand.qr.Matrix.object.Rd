\name{expand.qr.Matrix.object}
\alias{expand.qr.Matrix.object}
\title{
Expanded QR Decomposition Object
}
\description{
These are objects of class
"expand.qr.Matrix"
They represent the QR decomposition of a matrix.
}
\section{GENERATION}{
This class of objects is returned from the \code{expand.qr.Matrix} function
to represent a QR decomposition.
}
\section{METHODS}{
The following components must be included in a legitimate \code{"expand.qr.Matrix"}
object:
\item{q}{
The orthogonal factor.
}
\item{r}{
The factor containing an upper triangular factor portion.
}
\item{permutation}{
The column permutation applied to the original matrix in order to obtain the
factorization.
}}
\seealso{
\code{\link{qr.Matrix}}, \code{\link{expand}}
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
