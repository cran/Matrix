\name{Hermitian.test}
\alias{is.Hermitian}
\alias{Hermitian.test}
\title{
    Test for Symmetry or Conjugate Symmetry in a Matrix
}
\description{
    \code{Hermitian.test} returns the maximum modulus of the difference
    between a \code{Matrix} and its transpose.
    \code{is.Hermitian} checks if this value is less than a tolerance.
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
	tolerance for the maximum modulus of the elements of the difference 
	elements in the lower triangle and the corresponding elements in
	the upper triangle.
    }
}
\value{
    For \code{is.Hermitian}, a logical value is returned.
    If \code{tol} is missing, this value is equivalent to
    \code{inherits(x, "Hermitian")}.
    Otherwise the result is \code{TRUE} if the maximum modulus for the
    difference between elements in the lower triangle and the
    corresponding elements in the upper triangle is less than or equal
    to \code{tol}, and \code{FALSE} otherwise.
    For \code{Hermitian.test}, the maximum modulus for the difference
    between elements in the lower triangle and the corresponding
    elements in the upper triangle is returned.
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
