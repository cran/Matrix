\name{det.lu.Hermitian}
\alias{det.lu.Hermitian}
\title{
Determinant of a Hermitian Matrix from Triangular Factorization
}
\description{
Computes the determinant (or its logarithm) of a real symmetric or complex 
Hermitian Matrix given its symmetric indefinite factorization.
}
\usage{
det.lu.Hermitian(x, logarithm=T)
}
\arguments{
\item{x}{
An object of class \code{"lu.Hermitian"} representing the
symmetric indefinite decomposition of a real symmetric or
complex Hermitian matrix.
}
\item{logarithm}{
a logical variable indicating whether or not the logarithm of the modulus of
the determinant should be returned rather than the determinant itself.
The default is to return the logarithm.
}}
\value{
returns an object of class \code{"det"}.
}
\references{
Golub, G., and Van Loan, C. F. (1989).
\emph{Matrix Computations,}
2nd edition, Johns Hopkins, Baltimore.
}
\seealso{
\code{\link{det.object}}, \code{\link{det.Hermitian}}
}
\examples{
x <- Matrix( sample(-3:3, size = 9, replace = T), nrow = 3, ncol = 3)
x[row(x) > col(x)] <- t(x)[row(x) > col(x)]  # construct symmetric matrix
class(x) <-  Matrix.class(x)
det(lu(x))
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
