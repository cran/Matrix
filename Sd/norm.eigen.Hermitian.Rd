\name{norm.eigen.Hermitian}
\alias{norm.eigen.Hermitian}
\title{
Spectral Norm for Hermitian Matrices
}
\description{
Computes the spectral (2) norm or a real symmetric or complex Hermitian matrix
given its eigenvalue decomposition.
}
\usage{
norm.eigen.Hermitian(x)
}
\arguments{
\item{x}{
an object of class \code{"eigen.Hermitian"} representing
the eigenvalue decomposition of a real symmetric or complex Hermitian matrix.
}}
\value{
A numeric value of class \code{"norm"} 
representing the spectral or 2 norm of the matrix underlying \code{x}.
A copy of the call to \code{"norm"} is returned as an attribute.
}
\details{
The spectral or 2 norm the matrix underlying \code{x} is equal to its
eigenvalue of largest magnitude.
}
\references{
Golub, G., and Van Loan, C. F. (1989).
\emph{Matrix Computations,}
2nd edition, Johns Hopkins, Baltimore.
}
\seealso{
\code{\link{eigen}}, \code{\link{norm.Matrix}}, \code{\link{norm.svd.Matrix}}
}
\examples{
x <- Matrix( sample(-3:3, size = 9, replace = T), nrow = 3, ncol = 3)
x[row(x) > col(x)] <- t(x)[row(x) > col(x)]  # construct symmetric matrix
class(x) <- Matrix.class(x)
norm(eigen(x))
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
