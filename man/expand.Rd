\name{expand}
\title{Expand a Decomposition into Factors}
\usage{
expand(x, \dots)
}
\alias{expand}
\alias{expand.default}
\description{
  Expands decompositions stored in compact form into factors.
}
\arguments{
  \item{x}{a matrix decomposition.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\value{
  The expanded decomposition, typically a list of matrix factors.
}
\details{
  This is a generic function with special methods for different types
  of decompositions.
  Use \code{\link{methods}("expand")} to list all the methods for the
  \code{expand} generic.
}
\note{
  Factors for decompositions such as \code{lu} and \code{qr} can be
  stored in a compact form.
  The function \code{expand} allows all factors to be fully expanded.
}
\seealso{
  %    \code{\link{expand.lu.Matrix}}, \code{\link{expand.lu.Hermitian}},
  %    \code{\link{expand.qr.Matrix}}, \code{\link{expand.eigen.Matrix}},
  \code{\link{facmul}}.
}
\examples{
(x <- Matrix(round(rnorm(9),2), 3, 3))
(ex <- expand(lux <- lu(x)))
}
\keyword{algebra}
