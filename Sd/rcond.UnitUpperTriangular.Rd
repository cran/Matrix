\name{rcond.LowerTriangular}
\alias{rcond.LowerTriangular}
\alias{rcond.UpperTriangular}
\alias{rcond.UnitLowerTriangular}
\alias{rcond.UnitUpperTriangular}
\title{
Condition Estimation for Triangular Matrices
}
\description{
Computes an estimate of the reciprocal condition number of a triangular matrix.
}
\usage{
rcond.LowerTriangular(x, one.norm=T)
rcond.UpperTriangular(x, one.norm=T)
rcond.UnitLowerTriangular(x, one.norm=T)
rcond.UnitUpperTriangular(x, one.norm=T)
}
\arguments{
\item{x}{
A triangular Matrix, that is a matrix inheriting from class \code{"Matrix"} and one of \code{"LowerTriangular"}, \code{"UpperTriangular"}, \code{"UnitLowerTriangular"}, or \code{"UnitUpperTriangular"}.
}
\item{one.norm}{
A logical variable which indicates computation of the 1-norm condition 
estimate if \code{TRUE} and the infinity-norm condition estimate otherwise.
}}
\value{
A numeric value of class \code{"rcond"}, representing the reciprocal one or infinity
norm condition estimate.
A copy of the call to \code{"rcond"} is returned as an attribute.
}
\details{
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
Golub, G., and Van Loan, C. F. (1989).
\emph{Matrix Computations,}
2nd edition, Johns Hopkins, Baltimore.
}
\seealso{
\code{\link{rcond.Matrix}}
}
\examples{
x <- Matrix( sample(-3:3, size = 9, replace = T), nrow = 3, ncol = 3)
x[row(x) > col(x)] <- 0        # construct triangular matrix
class(x) <- Matrix.class(x)
rcond(x, one.norm = T)
rcond(x, one.norm = F)
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
