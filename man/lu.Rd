\name{lu}
\title{Triangular Decomposition of a Square Matrix}
\usage{
lu(x, \dots)
}
\alias{lu}
\description{
  Computes triangular decompositions of square matrices.
}
\arguments{
  \item{x}{a matrix. No missing values or IEEE special values are allowed.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\value{
  an object of class \code{"LU"}, see \code{\link{LU-class}}; this is
  a representation of a triangular decomposition of \code{x}.
}
\details{
  This is a generic function with special methods for different types
  of matrices.  Use \code{\link{methods}("lu")} to list all the methods
  for the \code{lu} generic.
}
\references{
  Golub, G., and Van Loan, C. F. (1989).
  \emph{Matrix Computations,}
  2nd edition, Johns Hopkins, Baltimore.
}
\seealso{
  \code{\link{LU-class}} and \code{\link{expand}};
  \code{\link{qr}}, \code{\link{chol}}.
}
\examples{
x <- Matrix(rnorm(9), 3, 3)
lu(x)
}
\keyword{array}
\keyword{algebra}
