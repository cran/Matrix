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
    \code{LowerTriangular.test} (\code{UpperTriangular.test}) returns
    the maximum modulus of the strict upper (lower) triangle of x.
    \code{is.LowerTriangular} and \code{is.UpperTriangular} check if
    this value is less than a tolerance.
}
\arguments{
    \item{x}{
	a numeric or complex matrix.
    }
    \item{tol}{
	tolerance for the maximum modulus of the elements in the upper
	(or lower) triangle of the matrix.
    }
}
\value{
    \code{LowerTriangular.test} (\code{UpperTriangular.test}) returns
    the maximum modulus of the elements in the strict upper (lower)
    triangle of the matrix \code{x}.
    
    \code{is.LowerTriangular} (\code{is.UpperTriangular}) returns a
    logical value according to whether the result of
    \code{LowerTriangular.test} (\code{UpperTriangular.test}) is less
    than or equal to \code{tol}.
%If \code{tol} is missing, this value is equivalent to 
%\code{inherits(x, "LowerTriangular")} or \code{inherits(x, "UpperTriangular")}. 
}
\examples{
x <- hilbert(9)
x[row(x) > col(x)] <- 0
is.UpperTriangular(x, tol = 0)
UpperTriangular.test(x)
}
\keyword{algebra}
