\name{solve.qr.Matrix}
\alias{solve.qr.Matrix}
\title{
Solve and Pseudo-Inverse with QR Decomposition
}
\description{
Given the QR decomposition of a matrix,
either solves a system of linear equations with that matrix as coefficient
matrix (in the least-squares sense if the system is underdetermined)
or else computes the pseudo-inverse of the matrix.
}
\usage{
solve.qr.Matrix(a, b, tol=0)
}
\arguments{
\item{a}{
An object of class \code{qr.Matrix}, representing the QR decomposition of a matrix.
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
transpose if complex) of \code{a} is to be used in the solve or inverse operation. 
The default is to use the untransposed matrix. 
The method does not allow a solve using the transpose of an 
underdetermined system.
}}
\value{
If \code{A} is the matrix whose QR decomposition is represented by \code{a},
an object of class \code{}Matrix"' is returned that is a least-squares solution 
\code{x} to the system of equations \code{A \%*\% x = b} if \code{A} has at least as many rows 
as columns and a basic solution to the system of equations if \code{A} has fewer
rows than columns. If \code{b} is not supplied, the pseudo-inverse of \code{A} is
returned.
Attributes include a copy of the call to \code{solve},
and the one norm reciprocal condition estimate if \code{tol} is nonnegative.
}
\details{
Can be used for matrices that are rank-deficient, i.e., 
whose rank is less than their minimum dimension.
}
\references{
Golub, G. H., and C. F. Van Loan (1989),
\emph{Matrix Computations,}
2nd edition, Johns Hopkins, Baltimore.
}
\seealso{
\code{\link{rcond.qr.Matrix}}, \code{\link{solve}}, \code{\link{solve.Matrix}}, \code{\link{solve.svd.Matrix}}.
}
\examples{
m <- sample(1:9, 1); n <- sample(1:9, 1)
a <- Matrix( sample(-9:9, m*n, replace = T), nrow = m, ncol = n)   
b <- rnorm(m)
z <- qr(a)                               # QR decomposition of a
if (m > n) t(a) \%*\% (a \%*\% solve(a,b) - b) else a \%*\% solve(a,b) - b
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
