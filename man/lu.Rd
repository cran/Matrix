\name{lu}
\title{Triangular Decomposition of a Square Matrix}
\usage{
lu(x, \dots)
}
\alias{lu}
\alias{lu.default}
\alias{lu.Matrix}
\description{
    Computes triangular decompositions of matrices.
}
\arguments{
    \item{x}{
	a matrix. No missing values or IEEE special values are allowed.
    }
    \item{...}{
	additional arguments, such as choosing lower or upper triangular factors
	in the Hermitian case.
    }
}
\value{
    a representation of a triangular decomposition of \code{x}.
}
\details{
    This is a generic function with special methods for different types
    of matrices.
    Use \code{\link{methods}("lu")} to list all the methods for the
    \code{lu} generic.
}

\section{BACKGROUND}{
    If \code{A} is a matrix, then \code{A = LU}, where \code{L} is unit
    lower trapezoidal, and \code{U} is upper trapezoidal.
    At least one of \code{L} or \code{U} is triangular; both are
    triangular if \code{A} is square. 

    If \code{A} is a Hermitian matrix, then \code{A} can be factored as
    \code{L B t(L)} or \code{t(L) B L}, where \code{L} is unit lower
    triangular, and \code{B} is block diagonal with either 1 by 1 or 2
    by 2 diagonal blocks.
    If \code{A} is positive semi-definite,
    then \code{B} is diagonal, and this factorization is a variant of the Choleski
    decomposition.
}
\references{
Golub, G., and Van Loan, C. F. (1989).
\emph{Matrix Computations,}
2nd edition, Johns Hopkins, Baltimore.
}
\seealso{\code{\link{qr}}, \code{\link{chol}}}
\examples{
x <- Matrix( rnorm(9), 3, 3)
\dontrun{
lu(x)
}
}
\keyword{array}
\keyword{algebra}
