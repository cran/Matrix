\name{det.qr.Matrix}
\alias{det.qr.Matrix}
\title{
Determinant of a Matrix from QR Decomposition
}
\description{
Computes the modulus of the determinant (or its logarithm) of a square Matrix  
from its QR decomposition.
}
\usage{
det.qr.Matrix(x, logarithm=T)
}
\arguments{
\item{x}{
an object of class \code{"qr.Matrix"} representing the QR decomposition of a square
Matrix.
}
\item{logarithm}{
a logical variable indicating whether or not the logarithm of the modulus of
the determinant should be returned rather than the determinant itself.
The default is to return the logarithm.
}}
\value{
returns an object of class \code{"det"}.
}
\details{
The determinant is computed as the the determinant of the triangular factor
in the QR factorization.
The sign of the determinant is not given since it involves computing
the determinant of the orthogonal or unitary factor, which may be less
efficient than computing the determinant from \code{x} directly.
}
\references{
Golub, G., and Van Loan, C. F. (1989).
\emph{Matrix Computations,}
2nd edition, Johns Hopkins, Baltimore.
}
\seealso{
\code{\link{det.object}}, \code{\link{det.Matrix}}.
}
\examples{
x <- Matrix( sample(-3:3, size = 9, replace = T), nrow = 3, ncol = 3)
det(qr(x))
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
