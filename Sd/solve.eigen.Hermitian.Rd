\name{solve.eigen.Hermitian}
\alias{solve.eigen.Hermitian}
\title{
Solve and Inverse with Hermitian Eigenvalue Decomposition
}
\description{
Given the eigenvalue decomposition of a real symmetric or complex Hermitian
matrix, either solves a system of linear equations
with that matrix as coefficient
matrix or else computes the inverse of the matrix.
}
\usage{
solve.eigen.Hermitian(a, b, tol=0)
}
\arguments{
\item{a}{
An object of class \code{eigen.Hermitian}, representing the symmetric 
indefinite decomposition of a real symmetric or complex Hermitian matrix.
}
\item{b}{
A real or complex matrix or vector. The number of rows in \code{b} must equal the
dimension of the matrix underlying \code{a}.
}
\item{tol}{
Tolerance for reciprocal condition number. Eigenvalues are considered
nonzero only if their ratio with the largest singular value exceeds \code{tol}.
By default, \code{tol} = 0.
}}
\value{
If \code{A} is the matrix whose eigenvalue decomposition is represented
by \code{a}, an object of class \code{"Matrix"} is returned that is the solution
\code{x} to the system of equations \code{A \%*\% x = b}
If \code{b} is not supplied, the inverse of \code{A} is returned.
The minimum least-squares solution or pseudo-inverse are returned if \code{A} is
rank deficient.
Attributes include a copy of the call to \code{solve}.
}
\details{
Can be used for matrices that are rank-deficient, i.e., whose rank
is less than their minimum dimension.
}
\references{
Golub, G. H., and C. F. Van Loan (1989),
\emph{Matrix Computations,}
2nd edition, Johns Hopkins, Baltimore.
}
\seealso{
\code{\link{rcond.eigen.Hermitian}}, \code{\link{solve}}, \code{\link{solve.Hermitian}}, \code{\link{solve.lu.Hermitian}}.
}
\examples{
 n <- 5
 a <- Matrix( rnorm(n*n), nrow = n, ncol = n)
 a[row(a) > col(a)] <- t(a)[row(a) > col(a)]  # construct symmetric matrix
 class(a) <- Matrix.class(x)
 b <- rnorm(n)
 z <- eigen(a)                                # eigenvalue decomposition
 a \%*\% solve(a,b) - b                         # residual
(solve(a) \%*\% b) - solve(a,b) 
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
