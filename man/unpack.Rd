\name{unpack}
\title{Full Storage Representation of Packed Matrices}
\usage{
unpack(x, \dots)
}
\alias{unpack}
\alias{unpack.default}
\description{
    Expands decompositions stored in compact form into matrix factors.
}
\arguments{
    \item{x}{
	a matrix stored in packed form. 
    }
    \item{...}{
	optional additional arguments (none are used in currently
	available methods)
    }
}
\value{
    A \code{Matrix} object containing the full-storage representation of \code{x}.
}
\details{
    This is a generic function with special methods for different types
    of packed matrices.
    Use \code{\link{methods}("unpack")} to list all the methods for the
    \code{unpack} generic.
}
\examples{
\dontrun{
x <- Diagonal( 1:3)
x
unpack(x)
}
}
\keyword{array}
\keyword{algebra}
