\name{facmul.lu.Hermitian}
\alias{facmul.lu.Hermitian}
\title{
Multiplication by Factors from a Symmetric Indefinite Decomposition
}
\description{
Performs multiplication by any of the factors in a symmetric or Hermitian
indefinite decomposition.
}
\usage{
facmul.lu.Hermitian(x, factor, y, transpose = F, left = T)
}
\arguments{
\item{x}{
an object of class \code{"lu.Hermitian"}.
}
\item{factor}{
either "T", "B", or "P",  where "T" selects the unit triangular factor
trapezoidal factor, "B" selects the block-diagonal factor, and "P selects the 
row permutation applied symmetrically to the original matrix in order to 
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
a logical value telling whether or not to multiply by the transpose of the
factor. The default is multiplication by the factor itself.
}
\item{left}{
a logical value telling whether the multiplication should occur on the left 
of \code{y} or on the right. The default is multiplication on the left.
}}
\value{
An object of class "Matrix" corresponding to the desired matrix product.
}
\details{
The triangular and block-diagonal factors for a symmetric indefinite 
decomposition from Lapack (Anderson et al., 1992) are stored in a 
single array the size of the underlying matrix. This functions allows these, 
as well as the row permutation of the original matrix that was used to form 
the factorization, to be applied separately without their explicit formation.
}
\references{
Anderson, E., et al. (1994).
\emph{LAPACK User's Guide,}
2nd edition, SIAM, Philadelphia.
}
\seealso{
\code{\link{lu.Hermitian}}, \code{\link{expand.lu.Hermitian}}
}
\examples{
 x <- Matrix( rnorm(81), nrow = 9, ncol = 9)
 x[row(x) > col(x)] <- t(x)[row(x) > col(x)]   # form symmetric matrix
 class(x) <- Matrix.class(x)
 z <- lu(x)                  # symmetric-indefinite factorization of x
 prod1 <- facmul(z,"T",facmul(z,"B", facmul(z,"T",transpose = T)))
 prod2 <- facmul(z,"P",facmul(z,"P",x,transpose = T,left=F))
 max(abs(prod1 - prod2))  # test product T B t(T) == P x t(P)
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
