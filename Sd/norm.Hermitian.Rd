\name{norm.Hermitian}
\alias{norm.Hermitian}
\title{
Norm of a Hermitian Matrix
}
\description{
Computes either the one or infinity norm, the Frobenius norm, or the maximum 
modulus among elements of a real symmetric or complex Hermitian matrix.
}
\usage{
norm.Hermitian(x, type="M", lower=T)
}
\arguments{
\item{x}{
a real symmetric or complex Hermitian Matrix inheriting from classes
\code{"Hermitian"} and \code{"Matrix"}.
}
\item{type}{
A character indication the type of norm desired. \code{"1"} or \code{"I"} indicates the one 
or infinity norm (maximum absolute column or row sum), \code{"F"} indicates the 
Frobenius seminorm (the Euclidean norm of \code{x} treated as if it were a vector),
and \code{"M"} indicates the maximum modulus of all the elements in \code{x}. 
The default is \code{"M"}
}
\item{lower}{
A logical variable telling whether to use the lower or upper triangle of the
matrix for the computation. The default is to select the lower triangle.
}}
\value{
A numeric value of class \code{"norm"}, 
representing the quantity chosen according to \code{type}.
A copy of the call to \code{"norm"} is returned as an attribute.
}
\details{
Based on the functions dlansy and zlanhe from Lapack (Anderson et al. (1994)).
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
\code{\link{Matrix}}, \code{\link{Matrix.class}}, \code{\link{norm}}. 
}
\examples{
x <- Matrix( sample(-3:3, size = 9, replace = T), nrow = 3, ncol = 3)
x[row(x) > col(x)] <- t(x)[row(x) > col(x)]  # construct symmetric matrix
class(x) <- Matrix.class(x)
norm(x, "1")
norm(x, "I")
norm(x, "F")
norm(x, "M")
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
