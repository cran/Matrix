\name{det.lu.Matrix}
\alias{det.lu.Matrix}
\title{
Determinant of a Matrix from LU Decomposition
}
\description{
Computes the determinant (or its logarithm) of a Matrix from its LU 
decomposition.
}
\usage{
det.lu.Matrix(x, logarithm=T)
}
\arguments{
\item{x}{
an object of class \code{"lu.Matrix"} representing the LU decomposition of a
square Matrix.
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
The determinant is computed from a triangular factorization obtained via the
function dgetrf or zgetrf from LAPACK (Anderson et al. 1994).
}
\references{
Anderson, E., et al. (1994).
\emph{LAPACK User's Guide,}
2nd edition, SIAM, Philadelphia.
Golub, G., and Van Loan, C. F. (1989).
\emph{Matrix Computations,}
2nd edition, Johns Hopkins, Baltimore.
}
\seealso{
\code{\link{det.object}}, \code{\link{det.Matrix}}.
}
\examples{
x <- Matrix( sample(-3:3, size = 9, replace = T), nrow = 3, ncol = 3)
det(lu(x))
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
