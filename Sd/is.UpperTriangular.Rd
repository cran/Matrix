\name{is.LowerTriangular}
\alias{is.LowerTriangular}
\alias{is.UpperTriangular}
\alias{LowerTriangular.test}
\alias{UpperTriangular.test}
\title{
Test for Triangularity in a Matrix
}
\description{
Either tests whether or not a Matrix is lower or upper triangular, or 
else gives a metric by which such a test can be made.
}
\usage{
is.LowerTriangular(x, tol)
is.UpperTriangular(x, tol)
LowerTriangular.test(x)
UpperTriangular.test(x)
}
\arguments{
\item{x}{
an numeric or complex matrix.
}
\item{tol}{
tolerance for the maximum modulus of the elements in the upper or lower
triangle of the matrix.
}}
\value{
For \code{is.LowerTriangular} and \code{is.UpperTriangular}, a logical value is returned.
If \code{tol} is missing, this value is equivalent to 
\code{inherits(x, "LowerTriangular")} or \code{inherits(x, "UpperTriangular")}. 
Otherwise the result is \code{TRUE} if the maximum modulus for the elements in the 
upper (repectively, lower) triangle is less than or equal to \code{tol} and 
\code{FALSE} otherwise.
For \code{LowerTriangular.test} and \code{UpperTriangular.test}, the maximum modulus for
the in the upper (respectively, lower) triangle is returned.
}
\examples{
x <- Matrix(1:9,3,3)
x[row(x) > col(x)] <- 0
is.UpperTriangular(x)
# [1] F
is.UpperTriangular(x, tol = 0)
# [1] T
UpperTriangular.test(x)
# [1] 0
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
