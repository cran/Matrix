\name{rcond.eigen.Hermitian}
\alias{rcond.eigen.Hermitian}
\title{
Condition Number for Hermitian Matrices from Eigenvalues
}
\description{
Computes the reciprocal two-norm condition number for a real symmetric
and complex Hermitian matrix from its eigenvalues.
}
\usage{
rcond.eigen.Hermitian(x)
}
\arguments{
\item{x}{
An object of class \code{"eigen.Hermitian"} representing
the eigenvalue decomposition of a real symmetric or complex Hermitian matrix.
}}
\value{
A numeric value of class \code{"rcond"}, representing the reciprocal two norm
condition number of the matrix underlying \code{x}.
A copy of the call to \code{"rcond"} is returned as an attribute.
}
\details{
The two-norm condition number of a symmetric or Hermitian matrix is the
ratio of smallest magnitude eigenvalue to its largest magnitude eigenvalue.
The condition number of a square matrix is the product of the norm of that
matrix and the norm of its inverse. Its values fall in the range \code{[1, Inf)},
where a value of \code{Inf} would imply a singular matrix. A matrix is said to
be ill-conditioned if its has a large condition number. Another way to
view a condition number of a matrix is as a factor by which errors for 
solutions to systems of equations with that matrix as coefficient matrix
can be multiplied. Condition numbers usually are estimated rather than 
computed exactly for reasons of efficiency. 
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
\code{\link{rcond.Hermitian}}
}
\examples{
x <- Matrix( sample(-3:3, size = 9, replace = T), nrow = 3, ncol = 3)
x[row(x) > col(x)] <- t(x)[row(x) > col(x)]  # construct symmetric matrix
class(x) <- Matrix.class(x)
rcond(eigen(x))
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
