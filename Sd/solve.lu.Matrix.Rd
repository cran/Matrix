\name{solve.lu.Matrix}
\alias{solve.lu.Matrix}
\title{
Solve and Inverse with LU Decomposition
}
\description{
Given the LU decomposition of a square matrix,
either solves a system of linear equations
with that matrix as coefficient 
matrix or else computes the inverse of the matrix.
}
\usage{
solve.lu.Matrix(a, b, tol=0, transpose=F, norma, tune, workspace)
}
\arguments{
\item{a}{
An object of class \code{lu.Matrix}, representing the LU decomposition of a
square matrix.
}
\item{b}{
A matrix or vector. If \code{transpose=T}, the number of rows of \code{b} must equal the
number of rows of the matrix underlying \code{a}, wheras if \code{transpose=F} the number
of rows of \code{b} must equal the number of columns of the matrix underlying \code{a}.
}
\item{tol}{
Tolerance for reciprocal condition estimation. If \code{tol} is negative, no
condition estimation is done. Otherwise, the reciprocal one norm condition
estimate is computed and the solve or inverse computation is done only if
the condition estimate is greater than \code{tol}. By default \code{tol = 0}.
}
\item{transpose}{
A logical variable indicating whether or not the transpose (conjugate
transpose if complex) of \code{a} is to be used in the solve or inverse operation. 
The default is to use the untransposed matrix.
}
\item{norma}{
The one norm of the matrix for use in condition estimation. By default,
it is assumed that the norm is available as an attribute of \code{a}, since the
default for \code{lu.Matrix} is to return the one and infinity norms of the matrix
in an attribute.
}
\item{tune}{
a integer vector or list of named tuning parameters that may affect
computational efficiency. The relevant parameters are the blocksize parameters
\code{NB}, \code{NBMIN}, and \code{NX} as described in \code{.laenv}.
}
\item{workspace}{
workspace provided to the underlying software.
The default is to use the optimum value relative to the tuning parameters.
The optimal workspace for the problem (with the given tuning parameters) is
included as part of the output attributes.
}}
\value{
If \code{A} is the matrix whose LU decomposition is represented 
by \code{a}, an object of class \code{}Matrix"' is returned that is the solution
\code{x} to the system of equations \code{A \%*\% x = b} 
If \code{b} is not supplied, the inverse of \code{A} is returned.
Attributes include a copy of the call to \code{solve},
the optimal workspace for the underlying software,
and the one norm reciprocal condition estimate if \code{tol} is nonnegative.
}
\details{
Based on the functions dgecon, dgetri, zgecon, and zgetri from
LAPACK (Anderson et al. 1994).


This function will not attempt a solve for \code{lu.Matrix} objects for 
matrices that are not square.
}
\references{
Anderson, E., et al. (1994).
\emph{LAPACK User's Guide,}
2nd edition, SIAM, Philadelphia.
}
\seealso{
\code{\link{solve}}, \code{\link{solve.lu.Hermitian}}, \code{\link{solve.Matrix}}, \code{\link{solve.qr.Matrix}}, \code{\link{solve.svd.Matrix}}.
}
\examples{
n <- 5
a <- Matrix( rnorm(n*n), nrow = n, ncol = n)z
b <- rnorm(n)
z <- lu(a)                                    # LU decomposition of a
a \%*\% solve(z,b) - b                          # residual
solve(z) \%*\% a                                # should be identity
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
