\name{LowerTriangular.test}
\title{Test a Matrix for Triangularity}
\usage{
is.LowerTriangular(x, tol = 0)
is.UpperTriangular(x, tol = 0)
LowerTriangular.test(x)
UpperTriangular.test(x)
}
\alias{is.LowerTriangular}
\alias{is.UpperTriangular}
\alias{LowerTriangular.test}
\alias{UpperTriangular.test}
\description{
    Either tests whether or not a Matrix is lower or upper triangular, or 
    else gives a metric by which such a test can be made.
}
\arguments{
    \item{x}{
	an numeric or complex matrix.
    }
    \item{tol}{
	tolerance for the maximum modulus of the elements in the upper or lower
	triangle of the matrix.
    }
}
\value{
    \code{LowerTriangular.test} (\code{UpperTriangular.test}) returns
    the maximum modulus of the elements in the strict upper (lower)
    triangle of the matrix \code{x}.
    
    For \code{is.LowerTriangular} (\code{is.UpperTriangular}) returns a
    logical value according to whether the result of
    \code{LowerTriangular.test} (\code{UpperTriangular.test}) is less
    than or equal to \code{tol}.
%If \code{tol} is missing, this value is equivalent to 
%\code{inherits(x, "LowerTriangular")} or \code{inherits(x, "UpperTriangular")}. 
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
