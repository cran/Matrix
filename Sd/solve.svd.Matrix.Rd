\name{solve.svd.Matrix}
\alias{solve.svd.Matrix}
\title{
Solve/Pseudo-Inverse with Singular Value Decomposition
}
\description{
Given the singular-value decomposition of a matrix,
either solves (in the least-squares sense) a system of linear equations 
with that matrix as coefficient
matrix or else computes the pseudo-inverse of the matrix.
}
\usage{
solve.svd.Matrix(a, b, tol=0, transpose=F)
}
\arguments{
\item{a}{
An object of class \code{svd.Matrix}, representing the singular-value decomposition
of a matrix.
}
\item{b}{
A matrix or vector. If \code{transpose=T} The number of rows of \code{b} must equal the
number of rows of the matrix underlying \code{a}, while if \code{transpose=F} the number
of rows of \code{b} must equal the number of columns of the matrix underlying \code{a}.
}
\item{tol}{
Tolerance for reciprocal condition number. Singular values are considered
nonzero only if their ratio with the largest singular value exceeds \code{tol}.
By default, \code{tol} = 0.
}
\item{transpose}{
A logical variable indicating whether or not the transpose (conjugate 
transpose if complex) of \code{a} is to be used for the solve or inverse operation. 
The default is to use the untransposed matrix.
}}
\value{
If \code{A} is the matrix whose singular-value decomposition is represented by \code{a},
an object of class \code{"Matrix"} is returned that is the minimum l2 norm
least-squares solution \code{x} to the system of equations \code{A \%*\% x = b} if \code{b} is
present; otherwise the pseudo-inverse of \code{A} is returned.
Attributes include a copy of the call to \code{solve}, the rank of the matrix
assumed in order to obtain the solution,
and the infinity norm reciprocal condition estimate if \code{tol} is nonnegative.
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
\code{\link{rcond.svd.Matrix}}, \code{\link{solve}}, \code{\link{solve.Matrix}}, \code{\link{solve.qr.Matrix}}.
}
\examples{
m <- sample(1:9, 1); n <- sample(1:9, 1)
a <- Matrix( sample(-9:9, m*n, replace = T), nrow = m, ncol = n)   
b <- rnorm(m)
z <- svd(a)                               # singular-value decomp of a
t(a) \%*\% (a \%*\% solve(a,b) - b)           # residual for normal equations
(solve(a) \%*\% b) - solve(a,b)
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
