\name{is.Orthonormal}
\alias{is.Orthonormal}
\alias{is.ColOrthonormal}
\alias{is.RowOrthonormal}
\alias{Orthonormal.test}
\title{
Test for Orthonormality in a Matrix
}
\description{
Either tests whether or not a Matrix has orthonormal rows and columns, or 
else gives a metric by which such a test can be made.
}
\usage{
is.Orthonormal(x, tol, byrow = F)
is.ColOrthonormal(x, tol)
is.RowOrthonormal(x, tol)
Orthonormal.test(x)
}
\arguments{
\item{x}{
an numeric or complex matrix.
}
\item{tol}{
tolerance for the testing for orthonormality.
}
\item{byrow}{
tells whether or not to do the tests using the rows rather than on the columns.
The default is to use the columns.
}}
\value{
For \code{is.Orthonormal}, \code{is.ColOrthonormal} and \code{is.RowOrthonormal} a 
logical value is returned.
If \code{tol} is missing from \code{is.Orthonormal}, this value is equivalent to 
\code{inherits(x, "Orthonormal")}. 
Otherwise the result is \code{TRUE} if the maximum of the modulus of the inner 
products between rows or columns of the matrix and the difference of the
norms of rows or columns from 1 is less than or equal to \code{tol},
and \code{FALSE} otherwise.


For \code{Orthonormal.test}, the maximum of the modulus of the
inner products between rows or columns of the matrix and the difference of the
norms of rows or columns from 1 is returned.
}
\seealso{
\code{\link{.Machine}}.
}
\examples{
v <- rnorm(3)
v <- v /vecnorm(v)
x <- diag(3) - 2 * (v \%*\% t(v))
is.Orthonormal(x)
# [1] F


is.Orthonormal(x, tol = .Machine$double.eps)
# [1] F


is.Orthonormal(x, tol = sqrt(.Machine$double.eps))
# [1] T


Orthonormal.test(x)
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
