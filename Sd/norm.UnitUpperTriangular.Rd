\name{norm.LowerTriangular}
\alias{norm.LowerTriangular}
\alias{norm.UpperTriangular}
\alias{norm.UnitLowerTriangular}
\alias{norm.UnitUpperTriangular}
\title{
Matrix Norm for Triangular Matrices.
}
\description{
Computes either the one or infinity norm, the Frobenius norm, 
or the maximum modulus among elements of a triangular matrix.
}
\usage{
norm.LowerTriangular(x, type = "M")
norm.UpperTriangular(x, type = "M")
norm.UnitLowerTriangular(x, type = "M")
norm.UnitUpperTriangular(x, type = "M")
}
\arguments{
\item{x}{
A triangular Matrix, that is, a Matrix inheriting from class
\code{"Matrix"} and one of \code{"LowerTriangular"}, \code{"UpperTriangular"},
\code{"UnitLowerTriangular"}, or \code{"UnitUpperTriangular"}.
}
\item{type}{
A character indication the type of norm desired. \code{"1"} indicates the one norm,
(maximum absolute column sum), \code{"I"} indicates the infinity norm (maximum
absolute row sum), \code{"F"} indicates the Frobenius norm (the Euclidean norm
of \code{x} treated as if it were a vector), and \code{"M"} indicates the maximum modulus
of all the elements in \code{x}. The default is \code{"M"}.
}}
\value{
A numeric value of class \code{"norm"}, representing the quantity chosen according 
to \code{type}.
A copy of the call to \code{"norm"} is returned as an attribute.
}
\details{
Based on the functions dlantr and zlantr from Lapack (Anderson et al. (1994)).
}
\references{
Anderson, E., et al. (1994).
\emph{LAPACK User's Guide,}
2nd edition, SIAM, Philadelphia.
Golub, G., and Van Loan, C. F. (1989).
\emph{Matrix Computations,}
2nd edition, Johns Hopkins, Baltimore.
}
\seealso{
\code{\link{norm.Matrix}}, \code{\link{vecnorm}}.
}
\examples{
x <- Matrix( sample(-3:3, size = 9, replace = T), nrow = 3, ncol = 3)
x[row(x) > col(x)] <- 0        # construct triangular matrix
class(x) <- Matrix.class(x)
norm(x, "1")
norm(x, "I")
norm(x, "F")
norm(x, "M")
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
