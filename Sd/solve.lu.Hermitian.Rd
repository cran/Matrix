\name{solve.lu.Hermitian}
\alias{solve.lu.Hermitian}
\title{
Solve and Inverse with Symmetric Indefinite Decomposition
}
\description{
Given the symmetric indefinite decomposition of a real symmetric 
or complex Hermitian matrix, either solves a system of linear equations
with that matrix as coefficient matrix or else computes the inverse of the 
matrix.
}
\usage{
solve.lu.Hermitian(a, b, tol=0, norma)
}
\arguments{
\item{a}{
An object of class \code{lu.Hermitian}, representing the 
symmetric indefinite decomposition of a real symmetric or
complex Hermitian matrix.
}
\item{b}{
A real or complex matrix or vector. The number of rows in \code{b} must equal the
dimension of the matrix underlying \code{a}.
}
\item{tol}{
Tolerance for reciprocal condition estimation. If \code{tol} is negative, no
condition estimation is done. Otherwise, the reciprocal one/infinity norm 
condition estimate is computed and the solve or inverse computation is done 
only if the condition estimate is greater than \code{tol}. By default, \code{tol = 0}.
}
\item{norma}{
The one/infinity norm of the matrix for use in condition estimation. By 
default, it is assumed that the norm is available as an attribute of \code{a},
since the default for \code{lu.Hermitian} is to return one/infinity norm of the
matrix as an attribute.
}}
\value{
If \code{A} is the matrix whose symmetric-indefinite decomposition is represented
by \code{a}, an object of class \code{"Matrix"} is returned that is the solution
\code{x} to the system of equations \code{A \%*\% x = b}
If \code{b} is not supplied, the inverse of \code{A} is returned.
Attributes include a copy of the call to \code{solve}, 
and the one/infinity norm reciprocal condition estimate if \code{tol} is 
nonnegative.
}
\details{
Based on the functions dsycon, dsytri, zhecon, zhetri from LAPACK 
(Anderson et al. 1994).
}
\references{
Anderson, E., et al. (1994).
\emph{LAPACK User's Guide,}
2nd edition, SIAM, Philadelphia.
}
\seealso{
\code{\link{rcond.Hermitian}}, \code{\link{solve}}, \code{\link{solve.Hermitian}}.
}
\examples{
 n <- 5
 a <- Matrix( rnorm(n*n), nrow = n, ncol = n)
 a[row(a) > col(a)] <- t(a)[row(a) > col(a)]  # construct symmetric matrix
 class(a) <- Matrix.class(x)
 b <- rnorm(n)
 z <- lu(a)                                   # symmetric-indefinite decomp
 a \%*\% solve(a,b) - b                         # residual
(solve(a) \%*\% b) - solve(a,b) 
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
