\name{det}
\title{Calculate the Determinant of a Matrix}
\alias{determinant}
\alias{determinant.Matrix}
\alias{determinant.UnitLowerTriangular}
\alias{determinant.LowerTriangular}
\alias{determinant.UnitUpperTriangular}
\alias{determinant.UpperTriangular}
\usage{
\method{determinant}{Matrix}(x, logarithm = TRUE, ...)
}
\description{
  \code{determinant} calculates the determinant (or its logarithm) of a matrix
}
\arguments{
  \item{x}{numeric matrix.}
  \item{logarithm}{logical.
    When \code{TRUE}, the default, the logarithm of the
    determinant is returned.
    When \code{FALSE}, the determinant of
    \code{x} is returned.
  }
  \item{\dots}{further arguments passed to or from other methods.}
}
\value{
 \code{determinant} return a list with two elements
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
  \code{determinant} is a generic function with special methods for different types
  of matrices or decompositions.
  Use \code{\link{methods}("determinant")} to list all the methods for the
  \code{determinant} generic.
  
}
\note{
  Often, computing the determinant is \emph{not} what you should be doing
  to solve a given problem.
}
\seealso{\code{\link[base]{det}}}
\examples{
h9 <- hilbert(9)
determinant(h9)
determinant(h9, log = FALSE)
}
\keyword{array}
\keyword{algebra}
