\name{rcond.lu.Hermitian}
\alias{rcond.lu.Hermitian}
\title{
Condition Estimation for Hermitian Matrices from Triangular Decomposition
}
\description{
Computes an estimate of the reciprocal condition number for a real symmetric
and complex Hermitian matrix from its symmetric indefinite decomposition.
}
\usage{
rcond.lu.Hermitian(x, normx)
}
\arguments{
\item{x}{
An object of class \code{"rcond.lu.Hermitian"} representing the symmetric 
indefinite decomposition real symmetric or complex Hermitian matrix.
}
\item{normx}{
The one/infinity norm of the matrix for use in condition estimation. By
default, it is assumed that the norm is available as an attribute of \code{a},
since the default for \code{lu.Hermitian} is to return one/infinity norm of the
matrix as an attribute.
}}
\value{
A numeric value of class \code{"rcond"}, representing the reciprocal one or infinity
norm condition estimate.
A copy of the call to \code{"rcond"} is returned as an attribute.
}
\details{
Based on the functions dsycon and zhecon from Lapack 
(Anderson et al. (1994)).
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
\code{\link{rcond.Matrix}}, \code{\link{rcond.Hermitian}}
}
\examples{
x <- Matrix( sample(-3:3, size = 9, replace = T), nrow = 3, ncol = 3)
x[row(x) > col(x)] <- t(x)[row(x) > col(x)]  # construct symmetric matrix
class(x) <- Matrix.class(x)
rcond(lu(x))
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
