\name{solve.LowerTriangular}
\alias{solve.LowerTriangular}
\alias{solve.UpperTriangular}
\alias{solve.UnitLowerTriangular}
\alias{solve.UnitUpperTriangular}
\title{
Solve and Inverse for Triangular Matrices
}
\description{
Either solves a system of linear equations with a triangular coefficient
matrix or else finds the inverse of a triangular matrix.
}
\usage{
solve.LowerTriangular(a, b, tol=0, transpose=F, left=T, scale.b=1)
solve.UpperTriangular(a, b, tol=0, transpose=F, left=T, scale.b=1)
solve.UnitLowerTriangular(a, b, tol=0, transpose=F, left=T, scale.b=1)
solve.UnitUpperTriangular(a, b, tol=0, transpose=F, left=T, scale.b=1)
}
\arguments{
\item{a}{
A triangular Matrix, inheriting from class \code{"Matrix"} and one of \code{"LowerTriangular"}, \code{"UpperTriangular"}, \code{"UnitLowerTriangular"}, or \code{"UnitUpperTriangular"}.
}
\item{b}{
A matrix or vector. The number of rows in \code{b} must equal the dimension of \code{a}.
}
\item{tol}{
Tolerance for reciprocal condition estimation. If \code{tol} is negative, no
condition estimation is done. Otherwise, the reciprocal infinity norm condition
estimate is computed and the solve or inverse computation is done only if
the condition estimate is greater than \code{tol}. By default, \code{tol = 0}.
}
\item{transpose}{
A logical variable indicating whether or not the transpose (conjugate 
transpose if complex) of \code{a} is to be used in the solve or inverse operation. 
The default is to use the untransposed matrix.
}
\item{left}{
A logical variable indicating whether or not the coefficient matrix \code{a} appears
on the left of the solution in the set of linear equations to be solved.
The default assumes that \code{a} appears to the left of the solution.
This argument is ignored when \code{b} is missing.
}
\item{scale.b}{
Scale factor to apply to the matrix \code{b}. The default is \code{scale.b = 1},
so that \code{b} is unscaled.
}}
\value{
An object of class \code{"Matrix"} that is the solution \code{x} to the system of 
equations \code{a \%*\% x = scale.b * b} if \code{b} is present and otherwise the inverse 
of \code{a}.
Attributes include a copy of the call to \code{solve},
and the one norm reciprocal condition estimate if \code{tol} is nonnegative.
}
\details{
Based on the functions dtrcon, dtrtrf, dtrtri, ztrcon, ztrtrf, and
ztrtri from LAPACK (Anderson et al. 1994).
}
\references{
Anderson, E., et al. (1994).
\emph{LAPACK User's Guide,}
2nd edition, SIAM, Philadelphia.
}
\seealso{
\code{\link{rcond.LowerTriangular}}, \code{\link{solve}}, \code{\link{solve.Matrix}}.
}
\examples{
 n <- 5
 a <- Matrix( rnorm(n*n), n, n)
 a[row(a) < col(a)] <- 0                 # form lower triangular matrix
 class(a) <- Matrix.class(x)
 b <- rnorm(n)
 a \%*\% solve(a,b) - b                    # residual
 solve(a) \%*\% a                          # should be identity
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
