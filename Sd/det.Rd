\name{det}
\alias{det}
\title{
    Determinant of a Matrix --- Generic Function
}
\description{
    Computes the determinant (or its logarithm) of a square matrix. 
    This is a generic function.
    Functions with names beginning in \code{det.} will be
    methods for this function.
    \code{"lu.Hermitian"}, \code{"lu.Matrix"}, \code{"qr.Matrix"}, \code{"schur.Matrix"},
    \code{"svd.Matrix"},  \code{"LowerTriangular"}, 
    \code{"UpperTriangular"}, \code{"UnitLowerTriangular"}, \code{"UnitUpperTriangular"}.
}
\usage{
det(x, logarithm = T, \dots)
}
\arguments{
    \item{x}{
	a square Matrix. No missing values or IEEE special values are allowed.
    }
    \item{logarithm}{
	logical value telling whether to compute the logarithm of the
	magnitude of the determinant, or the actual magnitude.
    }
    \item{...}{
	additional parameters, for example tuning constants for underlying code.
    }
}
\value{
    Returns an object of class \code{"det"}.
}
\section{BACKGROUND}{
    If \code{A} is a square matrix, then \code{det(A)} is the product of the eigenvalues of 
    \code{A}. It is not necessary to compute eigenvalues in order to compute a 
    determinant, however.
}
\references{
    Golub, G., and Van Loan, C. F. (1989).
    \emph{Matrix Computations,}
    2nd edition, Johns Hopkins, Baltimore.
}
\seealso{
    \code{\link{Matrix}}, \code{\link{det.Matrix}}, \code{\link{det.object}}.
}
\examples{
library(Matrix)
x <- Matrix( rnorm(9), 3, 3)
det(x) 
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.


