\name{det}
\title{Calculate the Determinant of a Matrix}
\alias{det}
\alias{det.default}
\alias{det.Matrix}
\alias{det.UnitLowerTriangular}
\alias{det.LowerTriangular}
\alias{det.UnitUpperTriangular}
\alias{det.UpperTriangular}
\usage{
det.Matrix(x, logarithm = TRUE)
}
\description{
    \code{det} calculates the determinant (or its logarithm) of a matrix
}
\arguments{
    \item{x}{numeric matrix.}
    \item{logarithm}{logical.
	When \code{TRUE}, the default, the logarithm of the
	determinant is returned.
	When \code{FALSE}, the determinant of
	\code{x} is returned.
    }
}
\value{
    \code{det.default} returns a numeric value.  Other \code{det}
    methods return a list with two elements
    \itemize{
	\item{modulus}{a numeric value.  The modulus (absolute value) of
	    the determinant or the logarithm of the modulus.  The value
	    of the \code{logarithm} argument is included as an
	    attribute.} 
	\item{sign}{a numeric value, which is \eqn{\pm 1}{+/- 1}
	    according to whether the determinant is positive or negative.}
    }
    }
\details{
    \code{det} is a generic function with special methods for different types
    of matrices or decompositions.
    Use \code{\link{methods}("det")} to list all the methods for the
    \code{det} generic.

    The \code{det} function from the base R package has a slightly
    different calling sequence.  In particular, it does not have a
    \code{logarithm} argument and always returns the determinant, never
    the logarithm of the determinant.  It is used as the default method for
    the generic \code{det} function defined here.  Be aware that the
    \code{logarithm} argument will only have an effect when \code{det}
    is called on objects that inherit from the \code{Matrix} class.
}
\note{
    Often, computing the determinant is \emph{not} what you should be doing
    to solve a given problem.
}
\examples{
h9 <- hilbert(9)
det(h9)
det(h9, log = FALSE)
}
\keyword{array}
\keyword{algebra}
