\name{det.eigen.Matrix}
\alias{det.eigen.Matrix}
\title{
Determinant of a Matrix from Eigenvalue Decomposition
}
\description{
Computes the determinant (or its logarithm) of a Matrix from its 
eigenvalue decomposition.
}
\usage{
det.eigen.Matrix(x, logarithm=T)
}
\arguments{
\item{x}{
an object of class \code{"eigen.Matrix"} representing
the eigenvalue decomposition of a Matrix.
}
\item{logarithm}{
a logical variable indicating whether or not the logarithm of the modulus of
the determinant should be returned rather than the determinant itself.
The default is to return the logarithm.
}}
\value{
returns an object of class \code{"det"}.
}
\references{
Golub, G., and Van Loan, C. F. (1989).
\emph{Matrix Computations,}
2nd edition, Johns Hopkins, Baltimore.
}
\seealso{
\code{\link{det.object}}, \code{\link{det.Hermitian}}, \code{\link{det.lu.Hermitian}}.
}
\examples{
x <- Matrix( sample(-3:3, size = 9, replace = T), nrow = 3, ncol = 3)
det(eigen(x))
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
