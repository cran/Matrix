\name{rcond}
\alias{rcond}
\title{
Reciprocal Condition Estimation for Matrices  --- Generic Function
}
\description{
Computes a variety of matrix norms and related functions (e.g., the element
of maximum modulus) for matrices. 
Classes which have methods for this function are available only in the Matrix
library and include: \code{"Matrix"}, \code{"LowerTriangular"}, and \code{"Hermitian"}.
}
\usage{
rcond(x, \dots)
}
\arguments{
\item{x}{
a matrix. No missing values or IEEE special values are allowed.
}
\item{...}{
additional arguments, for example to specify the matrix norm for which the
condition number is to be estimated.
}}
\value{
An estimate of the reciprocal condition number of \code{x}.
}
\section{BACKGROUND}{
The condition number of a matrix is the product of the matrix and the norm of 
its inverse (or pseudo-inverse if the matrix is not square). Since it can 
take on values between 1 and infinity, inclusive, it can be viewed as a 
measure of how close a matrix is to being rank deficient. It can also be
viewed as a factor by which errors in solving linear systems with this matrix
as coefficient matrix could be magnified. 


Condition numbers are usually estimated, since exact computation is costly
in terms of floating-point operations.
An (over) estimate of reciprocal condition number is given, since by doing
so overflow is avoided.
Matrices are well-conditioned if the reciprocal condition number is near 1
and ill-conditioned if it is near zero.
}
\references{
Golub, G., and Van Loan, C. F. (1989).
\emph{Matrix Computations,}
2nd edition, Johns Hopkins, Baltimore.
}
\seealso{
\code{\link{rcond.Matrix}}.
}
\examples{
library(Matrix)
x <- Matrix( rnorm(9), 3, 3)
rcond(x)                    
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
