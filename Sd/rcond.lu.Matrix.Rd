\name{rcond.lu.Matrix}
\alias{rcond.lu.Matrix}
\title{
Condition Estimation for Matrices from LU Decomposition
}
\description{
Computes an estimate of the reciprocal condition number of a matrix from its
LU decomposition.
}
\usage{
rcond.lu.Matrix(x, one.norm=T, normx)
}
\arguments{
\item{x}{
An object of class \code{"lu.Matrix"} representing the LU decomposition of a matrix.
}
\item{one.norm}{
A logical variable which indicateds computation of the 1-norm condition 
estimate if \code{TRUE} and the infinity-norm condition estimate otherwise.
}
\item{normx}{
The one or infinity norm of the matrix for use in computing the condition 
estimate. By default, it is assumed that the norm is available as an attribute
of \code{x}, since the default for \code{lu.Matrix} is to return the one and infinity 
norms of the matrix in an attribute.
}}
\value{
A numeric value of class \code{"rcond"}, representing the reciprocal one or infinity
norm condition estimate.
A copy of the call to \code{"rcond"} is returned as an attribute.
}
\details{
Based on the functions dgecon and zgecon from Lapack (Anderson et al. (1994)).
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
\code{\link{rcond.Matrix}}
}
\examples{
x <- Matrix( sample(-3:3, size = 9, replace = T), nrow = 3, ncol = 3)
rcond(lu(x), one.norm = T)
rcond(lu(x), one.norm = F)
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
