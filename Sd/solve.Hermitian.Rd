\name{solve.Hermitian}
\alias{solve.Hermitian}
\title{
Solve and Inverse for Hermitian Matrices
}
\description{
Either solves a system of linear equations with a real symmetric or
complex Hermitian coefficient matrix or else finds the inverse of a 
real symmetric or complex Hermitian matrix.
}
\usage{
solve.Hermitian(a, b, tol=0, lower=T, tune, workspace)
}
\arguments{
\item{a}{
A real symmetric or complex Hermitian Matrix, inheriting from class \code{"Hermitian"}.
}
\item{b}{
A real or complex matrix or vector. The number of rows in \code{b} must equal the
dimension of \code{a}.
}
\item{tol}{
Tolerance for reciprocal condition estimation. If \code{tol} is negative, no
condition estimation is done. Otherwise, the reciprocal one/infinity norm 
condition estimate is computed and the solve or inverse computation is done 
only if the condition estimate is greater than \code{tol}. By default, \code{tol = 0}.
}
\item{lower}{
A logical variable telling whether the lower or upper triangle of the matrix
is to be used for the computations. The default is to select the lower 
triangle.
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
An object of class \code{"Matrix"} that is the solution \code{x} to the system of
equations \code{a \%*\% x = scale.b * b} if \code{b} is present and otherwise the inverse
of \code{a}.
Attributes include a copy of the call to \code{solve}, an indicator telling which
triangle of the matrix was accessed for the computations, 
the optimal workspace for the underlying software,
and the one/infinity norm reciprocal condition estimate if \code{tol} is 
nonnegative.
}
\details{
Based on the functions dsycon, dsytrf, dsytri, zhecon, zhetrf, zhetri
from LAPACK (Anderson et al. 1994).
}
\references{
Anderson, E., et al. (1994).
\emph{LAPACK User's Guide,}
2nd edition, SIAM, Philadelphia.
}
\seealso{
\code{\link{rcond.Hermitian}}, \code{\link{solve}}, \code{\link{solve.lu.Hermitian}}, \code{\link{solve.eigen.Hermitian}}, \code{\link{solve.Matrix}}.
}
\examples{
 n <- 5
 a <- Matrix( rnorm(n*n), nrow = n, ncol = n)
 a[row(a) > col(a)] <- t(a)[row(a) > col(a)]  # construct symmetric matrix
 class(a) <- Matrix.class(x)
 b <- rnorm(n)
 a \%*\% solve(a,b) - b                         # residual
(solve(a) \%*\% b) - solve(a,b) 
}
\keyword{algebra}
\keyword{solve}
% Converted by Sd2Rd version 1.10.
