\name{facmul.lu.Matrix}
\alias{facmul.lu.Matrix}
\title{
Multiplication by Factors from an LU Decomposition
}
\description{
Performs multiplication by any of the factors in an LU decomposition.
}
\usage{
facmul.lu.Matrix(x, factor, y, transpose = F, left = T)
}
\arguments{
\item{x}{
an object of class \code{"lu.Matrix"}.
}
\item{factor}{
either "L", "U", or "P",  where "L" selects the unit lower-triangular or
trapezoidal factor, "U" selects the upper triangular or trapezoidal factor,
and "P selects the row permutation applied to the original matrix in order to 
obtain the decomposition.
}
\item{y}{
a numeric or complex vector or matrix. The default is the identity 
matrix of the order equal to the number of rows (\code{left = T}) or 
number of columns (\code{left = F}) in the matrix underlying \code{x}, so that either 
the factor itself (\code{transpose = F}) or its transpose 
(\code{transpose = T}) are returned when \code{y} is missing.
}
\item{transpose}{
a logical value telling whether or not to multiply by the (if complex,
conjugate) transpose of the factor. The default is multiplication by the
factor itself.
}
\item{left}{
a logical value telling whether the multiplication should occur on the left 
of \code{y} or on the right. The default is multiplication on the left.
}}
\value{
An object of class "Matrix" corresponding to the desired matrix product.
}
\details{
The lower and upper triangular factors for an LU decomposition from Lapack 
(Anderson et al., 1994) are stored in a single array the size of the 
underlying matrix. This functions allows these, as well as the permutation,
to be applied separately without their explicit formation.
}
\references{
Anderson, E., et al. (1994).
\emph{LAPACK User's Guide,}
2nd edition, SIAM, Philadelphia.
}
\seealso{
\code{\link{lu.Matrix}}, \code{\link{expand.lu.Matrix}}
}
\examples{
 x <- Matrix( rnorm(9), nrow = 3, ncol = 3)
 z <- lu(x)                                        # LU factorization of x
 max(abs(facmul(z,"L",facmul(z,"U")) - facmul(z,"P",x))) # test product LU
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
