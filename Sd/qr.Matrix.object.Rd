\name{qr.Matrix.object}
\alias{qr.Matrix.object}
\title{
QR Decomposition Object
}
\description{
These are objects of class
"qr.Matrix"
They represent the QR decomposition of a matrix.
}
\section{GENERATION}{
This class of objects is returned from the \code{qr.Matrix} function
to represent a QR decomposition.
}
\section{METHODS}{
The \code{"qr.Matrix"} class of objects has methods for the following generic
functions:
\code{det}, \code{expand}, \code{facmul}, \code{rcond}, \code{solve}
}
\section{STRUCTURE}{
The following component must be included in a legitimate \code{"qr.Matrix"} object:
\item{factors}{
a list composed of two unnamed components. The first element of this list 
is a matrix which has the R factor of the QR decomposition in its upper 
triangle and information on the transformations whose product make up the 
orthogonal Q factor in its lower triangle. The second element of the list is
a vector of auxiliary information needed to apply or reconstruct the
orthogonal or Q factor of the QR decomposition.
The following component may also be included in a \code{"qr.Matrix"} object:
}
\item{pivot}{
if pivoting is requested, this is the column permutation applied to the 
original matrix to achieve the QR decomposition.
}}
\seealso{
\code{\link{qr.Matrix}}
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
