\name{is.Hermitian}
\alias{is.Hermitian}
\alias{Hermitian.test}
\title{
Test for Symmetry or Conjugate Symmetry in a Matrix
}
\description{
Either tests whether or not a Matrix is real symmetric or complex Hermitian, 
or else gives a metric by which such a test can be made.
}
\usage{
is.Hermitian(x, tol)
Hermitian.test(x)
}
\arguments{
\item{x}{
an numeric or complex matrix.
}
\item{tol}{
tolerance for the maximum modulus for the difference between elements in the
lower triangle and the corresponding elements in the upper triangle.
}}
\value{
For \code{is.Hermitian}, a logical value is returned.
If \code{tol} is missing, this value is equivalent to \code{inherits(x, "Hermitian")}.
Otherwise the result is \code{TRUE} if the maximum modulus for the difference 
between elements in the lower triangle and the corresponding elements in the 
upper triangle is less than or equal to \code{tol}, and \code{FALSE} otherwise.
For \code{Hermitian.test}, the maximum modulus for the difference between elements 
in the lower triangle and the corresponding elements in the upper triangle is 
returned.
}
\examples{
x <- Matrix(1:9,3,3)
x[row(x) > col(x)] <- t(x)[row(x) > col(x)]
is.Hermitian(x)
# [1] F


is.Hermitian(x, tol = 0)
# [1] T


Hermitian.test(x)
# [1] 0
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
