\name{norm.svd.Matrix}
\alias{norm.svd.Matrix}
\title{
    Spectral Norm from Eigenvalue Decomposition
}
\description{
    Computes the spectral (2) norm of a matrix given its singular-value 
    decomposition.
}
\usage{
norm.svd.Matrix(x)
}
\arguments{
    \item{x}{
	An object of class \code{"svd.Matrix"} representing the 
	singular value decomposition of a numeric or complex Matrix.
    }
}
\value{
    An object of class \code{"norm"} representing the spectral or 2 norm
    of the matrix underlying \code{x}.
    A copy of the call to \code{norm} is returned as an attribute.
}
\details{
    The spectral or 2 norm of a matrix is equal to its largest singular value.
}
\references{
    Golub, G., and Van Loan, C. F. (1989).
    \emph{Matrix Computations,}
    2nd edition, Johns Hopkins, Baltimore.
}
\seealso{
    \code{\link{svd}}, \code{\link{norm.Matrix}},
    \code{\link{norm.eigen.Hermitian}}.
}
\examples{
x <- Matrix( sample(-3:3, size = 9, replace = T), nrow = 3, ncol = 3)
norm(svd(x))
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.








