\name{det.LowerTriangular}
\alias{det.LowerTriangular}
\alias{det.UpperTriangular}
\alias{det.UnitLowerTriangular}
\alias{det.UnitUpperTriangular}
\title{
Determinant of a Triangular Matrix
}
\description{
Computes the determinant (or its logarithm) of a triangular Matrix.
}
\usage{
det.LowerTriangular(x, logarithm=T)
det.UpperTriangular(x, logarithm=T)
det.UnitLowerTriangular(x, logarithm=T)
det.UnitUpperTriangular(x, logarithm=T)
}
\arguments{
    \item{x}{
	a triangular Matrix, that is, a Matrix inheriting from class
	\code{"Matrix"}, and from one of \code{"LowerTriangular"},
	\code{"UpperTriangular"}, \code{"UnitLowerTriangular"}, and
	\code{"UnitUpperTriangular"}.
    }
    \item{logarithm}{
	a logical variable indicating whether or not the logarithm of
	the modulus of the determinant should be returned rather than
	the determinant itself. 
	The default is to return the logarithm.
    }
}
\value{
    returns an object of class \code{"det"}.
}
\references{
    Golub, G., and Van Loan, C. F. (1989).
    \emph{Matrix Computations,}
    2nd edition, Johns Hopkins, Baltimore.
}
\seealso{
    \code{\link{det.object}}, \code{\link{det.Matrix}}.
}
\examples{
x <- Matrix( sample(-3:3, size = 9, replace = T), nrow = 3, ncol = 3)
x[row(x) > col(x)] <- 0        # construct triangular matrix
class(x) <-  Matrix.class(x)
det(x)
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
