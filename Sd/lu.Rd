\name{lu}
\alias{lu}
\title{
Triangular Decomposition of a Matrix --- Generic Function
}
\description{
Computes triangular decompositions of matrices.
This is a generic function.
Functions with names beginning in \code{lu.} will be
methods for this function.
\code{"Matrix"}, \code{"Hermitian"}
}
\usage{
lu(x, \dots)
}
\arguments{
\item{x}{
a matrix. No missing values or IEEE special values are allowed.
}
\item{...}{
additional arguments, such as choosing lower or upper triangular factors
in the Hermitian case.
}}
\value{
a representation of a triangular decomposition of \code{x}.
}
\section{BACKGROUND}{
If \code{A} is a matrix, then \code{A = LU}, where \code{L} is unit lower trapezoidal, and \code{U}
is upper trapezoidal. At least one of \code{L} or \code{U} is triangular; both are
triangular if `A" is square. 


If \code{A} is a Hermitian matrix, then \code{A} can be factored as \code{L B t(L)} or 
\code{t(L) B L}, where \code{L} is unit lower triangular, and \code{B} is block diagonal
with either 1 by 1 or 2 by 2 diagonal blocks. If \code{A} is positive semi-definite,
then \code{B} is diagonal, and this factorization is a variant of the Choleski
decomposition.
}
\references{
Golub, G., and Van Loan, C. F. (1989).
\emph{Matrix Computations,}
2nd edition, Johns Hopkins, Baltimore.
}
\seealso{
\code{\link{qr}}, \code{\link{chol}}
}
\examples{
library(Matrix)
x <- Matrix( rnorm(9), 3, 3)
lu(x) 
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
