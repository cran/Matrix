\name{Hermitian.test}
\title{Test a Matrix for Conjugate Symmetry}
\usage{
is.Hermitian(x, tol)
Hermitian.test(x)
}
\alias{is.Hermitian}
\alias{Hermitian.test}
\description{
    \code{Hermitian.test} returns the maximum modulus of the difference
    between a matrix and the conjugate of its transpose.
    \code{is.Hermitian} checks if this value is less than a tolerance.
}
\arguments{
    \item{x}{
	a numeric or complex matrix.
    }
    \item{tol}{
	tolerance for the maximum modulus of the elements of the difference 
	elements in the lower triangle and the corresponding elements in
	the upper triangle.
    }
}
\value{
    \code{Hermitian.test} returns the maximum modulus of the difference
    between \code{x} and its conjugate transpose.
    (If \code{x} is not a matrix or is not square, \code{Inf} is returned.)
    \code{is.Hermitian} returns a logical value according to whether
    the value of \code{Hermitian.test(x)} is less than or equal to the
    tolerance \code{tol}.
}
\examples{
is.Hermitian(hilbert(6))  # returns TRUE
x <- Matrix(1:9,3,3)
x[row(x) > col(x)] <- t(x)[row(x) > col(x)]
is.Hermitian(x, tol = 0)  # returns TRUE
Hermitian.test(x)         # returns 0
}
\keyword{algebra}
