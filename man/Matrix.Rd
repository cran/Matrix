\name{Matrix}
\alias{Matrix}
\title{Construct a Classed Matrix}
\usage{
Matrix(data=NA, nrow=1, ncol=1, byrow=FALSE, dimnames=NULL)
}
\description{
  Construct a Matrix of a class that inherits from \code{Matrix}.
}
\arguments{
  \item{data}{an optional numeric data vector or matrix.}
  \item{nrow}{the desired number of rows}
  \item{ncol}{the desired number of columns}
  \item{byrow}{logical.  If \code{FALSE} (the default) the matrix is
    filled by columns, otherwise the matrix is filled by rows.}
  \item{dimnames}{a \code{\link{dimnames}} attribute for the matrix: a
    \code{list} of two character components.}
}
\value{
  Returns an \code{nrow} by \code{ncol} matrix of a class that
  inherits from \code{"Matrix"}.
}
\details{
  If either of \code{nrow} or \code{ncol} is not given, an attempt is
  made to infer it from the length of \code{data} and the other
  parameter.

  Although it is sometime possible to mix unclassed matrices (created
  with \code{matrix}) with ones of class \code{"Matrix"}, it is much
  safer to always use carefully constructed ones of class
  \code{"Matrix"}.
}
\seealso{
  \code{\link{Matrix-class}}; \code{\link{matrix}}.
}
\examples{
Matrix(0, 3, 2)             # a 3 by 2 matrix of zeros
Matrix(1:6, 3, 2)           # a 3 by 2 matrix
Matrix(1:6, nrow=3)
Matrix(1:6, ncol=2)
Matrix(1:9, nrow=3,
       dimnames = list(c("a", "b", "c"), c("A", "B", "C")))
}
\keyword{array}
\keyword{algebra}
