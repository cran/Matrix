\name{det}
\title{Calculate the Determinant of a Matrix}
\alias{det}
\alias{det.Matrix}
\alias{det.UnitLowerTriangular}
\alias{det.LowerTriangular}
\alias{det.UnitUpperTriangular}
\alias{det.UpperTriangular}
\usage{
det(x, method = c("qr","eigenvalues"))
}
\description{
    \code{det} calculates the determinant (or its logarithm) of a matrix
}
\arguments{
    \item{x}{numeric matrix.}
    \item{log}{logical. When \code{FALSE}, the default, the determinant of
	\code{x} is returned.  When \code{TRUE}, the logarithm of the
	determinant is returned.
    }
}
\value{
    A numeric value with class \code{"det"}.
}
\note{
    Often, computing the determinant is \emph{not} what you should be doing
    to solve a given problem.
}
\examples{
x <- matrix(1:4, ncol=2)
\dontrun{
det(x)
det(x, log = TRUE)

det(print(cbind(1,1:3,c(2,0,1))))
}
}
\keyword{array}
\keyword{algebra}
