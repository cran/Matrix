\name{solve.Matrix}
\alias{solve.Matrix}
\title{
Solve and Inverse for General Matrices
}
\description{
Either solves a system of linear equations or else
finds the inverse or pseudo-inverse of a matrix.
}
\usage{
solve.Matrix(a, b, tol=0, transpose=F, tune, workspace)
}
\arguments{
\item{a}{
A numeric or complex Matrix, inheriting from class \code{"Matrix"}.
}
\item{b}{
A matrix or vector. If \code{transpose=T} The number of rows of \code{b} must equal the
number of rows of \code{a}, while if \code{transpose=F} the number of rows of \code{b} must
equal the number of columns of \code{a}.
}
\item{tol}{
Tolerance for reciprocal condition estimation. If \code{tol} is negative, no
condition estimation is done. Otherwise, the reciprocal one norm condition
estimate is computed and the solve or inverse computation is done only if
the condition estimate is greater than \code{tol}. By default \code{tol = 0}.
}
\item{transpose}{
A logical variable indicating whether or not the transpose (conjugate
transpose if complex) of \code{a} is to be used in the solve or inverse 
operation. The default is to use the untransposed matrix.
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
An object of class \code{}Matrix"'.
If \code{a} is square, it is the solution \code{x} to the system of
equations \code{a \%*\% x = scale.b * b} if \code{b} is present and otherwise the inverse
of \code{a}. 
If \code{a}is not square, it is the the minimum l-2 norm solution to the system 
of equations if \code{b} is present and otherwise the pseudo-inverse of \code{a}.
Attributes include a copy of the call to \code{solve},
the optimal workspace for the underlying software,
and the one norm reciprocal condition estimate if \code{tol} is nonnegative.
}
\details{
Based on the functions dgecon, dgetrf, dgetri, zgecon, zgetrf, zgetri,
dtrcon, dgeqpf, dormqr, dorgqr, ztrcon, zgeqpf, zunmqr, and zungqr from 
LAPACK (Anderson et al. 1994).
}
\references{
Anderson, E., et al. (1994).
\emph{LAPACK User's Guide,}
2nd edition, SIAM, Philadelphia.
}
\seealso{
\code{\link{rcond.Matrix}}, \code{\link{solve}}, \code{\link{solve.Hermitian}}, \code{\link{solve.lu.Matrix}}, \code{\link{solve.qr.Matrix}}, \code{\link{solve.svd.Matrix}}.
}
\examples{
m <- sample(1:9,1); n <- sample(1:9,1)
a <- Matrix( rnorm(m*n), nrow = m, ncol = n)
b <- rnorm(m)
t(a) \%*\% (a\%*\% solve(a,b) - b)            # residual of normal equations
(solve(a) \%*\% b) - solve(a,b)               
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
