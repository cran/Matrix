\name{rcond.qr.Matrix}
\alias{rcond.qr.Matrix}
\title{
Condition Estimation for Matrices from QR Decomposition
}
\description{
Computes an estimate of the reciprocal condition number of a matrix from its
QR decomposition.
}
\usage{
rcond.qr.Matrix(x, one.norm=T)
}
\arguments{
\item{x}{
An object of class \code{"qr.Matrix"}, representing the QR decomposition of a 
matrix.
}
\item{one.norm}{
A logical variable which indicateds computation of the 1-norm condition 
estimate if \code{TRUE} and the infinity-norm condition estimate otherwise.
}}
\value{
A numeric value of class \code{"rcond"}, representing the reciprocal one or infinity
norm condition estimate.
A copy of the call to \code{"rcond"} is returned as an attribute.
}
\details{
The condition estimate returned is that of the triangular factor, using 
complete orthogonalization if the underlying matrix has more columns than rows.
Based on the functions dtzrqf, dtrcon, ztzrqf, ztrcon from Lapack 
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
\code{\link{rcond.qr.Matrix}}, \code{\link{rcond.svd.Matrix}}
}
\examples{
x <- Matrix( sample(-3:3, size = 12, replace = T), nrow = 4, ncol = 3)
rcond(qr(x), one.norm = T)
rcond(qr(x), one.norm = F)
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
