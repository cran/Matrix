\name{facmul}
\alias{facmul}
\title{
Multiplication by Decomposition Factors --- Generic Function
}
\description{
Performs multiplication by factors for certain decompositions (and allows
explicit formation of those factors).
This is a generic function.
Functions with names beginning in \code{facmul.} will be
methods for this function.
\code{"lu.Matrix"}, \code{"lu.Hermitian"}, \code{"qr.Matrix"}.
}
\usage{
facmul(x, factor, y, transpose = F, left = T, \dots)
}
\arguments{
\item{x}{
a matrix decomposition. 
No missing values or IEEE special values are allowed.
}
\item{factor}{
an indicator for selecting a particular factor for multiplication.
}
\item{y}{
a matrix or vector to be multiplied by the factor or its transpose.
No missing values or IEEE special values are allowed.
}
\item{transpose}{
a logical value indicating whether or not the transpose of the factor 
should be applied. The default is to use the untransposed factor.
}
\item{left}{
a logical value indicating whether the multiplication should occur on
the left or on the right of \code{y}. The default is multiplication on the left.
}
\item{...}{
the method for \code{"qr.Matrix"} has additional arguments.
}}
\value{
the product of the selected factor (or its transpose) and \code{y}
}
\section{NOTE}{
Factors for decompositions such as \code{lu} and \code{qr} can be stored in a compact
form. The function \code{facmul} allows multiplication without explicit formation
of the factors, saving both storage and operations.
}
\references{
Golub, G., and Van Loan, C. F. (1989).
\emph{Matrix Computations,}
2nd edition, Johns Hopkins, Baltimore.
}
\seealso{
\code{\link{facmul.lu.Matrix}}, \code{\link{facmul.lu.Hermitian}}, \code{\link{facmul.qr.Matrix}}, \code{\link{expand}}.
}
\examples{
library(Matrix)
x <- Matrix(rnorm(9), 3, 3)
qrx <- qr(x)                      # QR factorization of x
y <- rnorm(3)
facmul( qr(x), factor = "Q", y)   # form Q y
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
