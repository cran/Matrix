\name{diagDet}
\alias{diagDet}
\title{Determinant of triangular matrices}
\description{
    Return the determinant (or log-determinant) of triangular matrices.
}
\usage{
diagDet(x, logarithm=TRUE)
}
\arguments{
    \item{x}{a numeric vector representing the diagonal of a triangular or
	a diagonal matrix.}
    \item{logarithm}{logical.
	When \code{TRUE}, the default, the logarithm of the
	determinant is returned.
	When \code{FALSE}, the determinant of
	\code{x} is returned.
    }
}
\details{
    The determinant of a triangular matrix can be calculated from its
    diagonal elements only.
}
\value{
    An object of class \code{det} as a list with two elements
    \itemize{
	\item{modulus}{a numeric value.  The modulus (absolute value) of
	    the determinant or the logarithm of the modulus.  The value
	    of the \code{logarithm} argument is included as an
	    attribute.} 
	\item{sign}{a numeric value, which is \eqn{\pm 1}{+/- 1}
	    according to whether the determinant is positive or negative.}
    }
}
\note{
    Often, computing the determinant is \emph{not} what you should be doing
    to solve a given problem.
}
\seealso{\code{\link{det}}}
\examples{
diagDet(1:5)  # log of product
diagDet(1:5, logarithm = FALSE)
}
\keyword{array}
\keyword{algebra}
