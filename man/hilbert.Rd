\name{hilbert}
\alias{hilbert}
\title{Generate a Hilbert matrix}
\description{
    Generate the \code{n} by \code{n} symmetric Hilbert matrix.  Because
    these matrices are ill-conditioned for moderate to large \code{n},
    they are often used for testing numerical linear algebra code.
}
\usage{
hilbert(n)
}
\arguments{
  \item{n}{a positive integer.}
}
\value{
    the \code{n} by \code{n} symmetric Hilbert matrix.
}
\examples{
round(hilbert(6), 4)
}
\keyword{array}
\keyword{algebra}
