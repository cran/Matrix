\name{det.Matrix}
\alias{det.Matrix}
\title{
Determinant of a Matrix
}
\description{
Computes the determinant (or its logarithm) of a Matrix.
}
\usage{
det.Matrix(x, logarithm=T, tune)
}
\arguments{
\item{x}{
a numeric or complex square Matrix inheriting from class \code{"Matrix"}.
}
\item{logarithm}{
a logical variable indicating whether or not the logarithm of the modulus of
the determinant should be returned rather than the determinant itself.
The default is to return the logarithm.
}
\item{tune}{
a integer vector or list of named tuning parameters that may affect
computational efficiency. The relevant parameter is the blocksize parameters
\code{NB} as described in \code{.laenv}.
}}
\value{
returns an object of class \code{"det"}.
}
\details{
The determinant is computed from an LU factorization obtained via the
function dgetrf or zgetrf from LAPACK (Anderson et al. 1994).
}
\references{
Anderson, E., et al. (1994).
\emph{LAPACK User's Guide,}
2nd edition, SIAM, Philadelphia.
Golub, G., and Van Loan, C. F. (1989).
\emph{Matrix Computations,}
2nd edition, Johns Hopkins, Baltimore.
}
\seealso{
\code{\link{Matrix}}, \code{\link{det.object}}, \code{\link{det.lu.Matrix}}, \code{\link{det.svd.Matrix}}.
}
\examples{
x <- Matrix( sample(-3:3, size = 9, replace = T), nrow = 3, ncol = 3)
det(x)
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
