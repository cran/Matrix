\name{norm.Matrix}
\alias{norm.Matrix}
\title{
    Norm of a Matrix
}
\description{
    Computes either the one or infinity norm, the Frobenius norm, 
    or the maximum modulus among elements of a matrix.
}
\usage{
norm.Matrix(x, type = "M")
}
\arguments{
    \item{x}{
	a real or complex Matrix inheriting from class \code{"Matrix"}.
    }
    \item{type}{
	A character indication the type of norm desired.
	\code{"1"} indicates the one norm, (maximum absolute column sum),
	\code{"I"} indicates the infinity norm (maximum absolute row sum),
	\code{"F"} indicates the Frobenius norm (the Euclidean norm
	of \code{x} treated as if it were a vector), and \code{"M"}
	indicates the maximum modulus of all the elements in
	\code{x}. The default is \code{"M"}.
    }
}
\value{
    A numeric value of class \code{"norm"}, representing the quantity
    chosen according to \code{type}.
    A copy of the call to \code{"norm"} is returned as an attribute.
}
\details{
    Based on the functions dlange and zlange from Lapack (Anderson et
    al. (1994)).
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
    \code{\link{norm.eigen.Matrix}}, \code{\link{norm.Hermitian}},
    \code{\link{norm.LowerTriangular}}, 
    \code{\link{norm.svd.Matrix}}.
}
\examples{
x <- Matrix(sample(-3:3, size = 9, replace = T), nrow = 3, 
            ncol = 3)
norm(x, "1")
norm(x, "I")
norm(x, "F")
norm(x, "M")
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
