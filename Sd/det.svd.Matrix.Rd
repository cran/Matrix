\name{det.svd.Matrix}
\alias{det.svd.Matrix}
\title{
Determinant of a Matrix from Singular-Value Decomposition
}
\description{
Computes the modulus of the determinant (or its logarithm) of a square Matrix 
from its singular-value decomposition.
}
\usage{
det.svd.Matrix(x, logarithm=T)
}
\arguments{
\item{x}{
an object of class \code{"svd.Matrix"} representing the singular-value 
decomposition of a square Matrix.
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
The sign of the determinant is returned as \code{NA} since it cannot be efficiently
computed from the singular value decomposition.
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
det(svd(x))
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
