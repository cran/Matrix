\name{det.Hermitian}
\alias{det.Hermitian}
\title{
Determinant of a Hermitian Matrix
}
\description{
Computes the determinant (or its logarithm) of a real symmetric or 
complex Hermitian Matrix.
}
\usage{
det.Hermitian(x, logarithm=T, lower=T, tune, workspace)
}
\arguments{
\item{x}{
a real symmetric or complex Matrix inheriting from the classes
\code{"Hermitian"} and \code{"Matrix"}.  
}
\item{logarithm}{
a logical variable indicating whether or not the logarithm of the modulus of
the determinant should be returned rather than the determinant itself.
The default is to return the logarithm.
}
\item{lower}{
a logical variable telling whether to use the lower or upper triangle of the
matrix for the factorization used to compute the determinant. The default is 
to use the lower triangle.
}
\item{tune}{
a integer vector or list of named tuning parameters that may affect
computational efficiency. The relevant parameters are the blocksize parameters
\code{NB} and \code{NBMIN} as described in \code{.laenv}.
}
\item{workspace}{
workspace provided to the underlying software.
The default is to use the optimum value relative to the tuning parameters.
The optimal workspace for the problem (for the tuning parameters) is
included as part of the output attributes.
}}
\value{
returns an object of class \code{"det"}.
}
\details{
The determinant is computed from a triangular factorization obtained via the 
functions dsytrf or zhetrf from LAPACK (Anderson et al. 1994).
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
\code{\link{det.object}}, \code{\link{det.Matrix}}, \code{\link{det.lu.Hermitian}}, \code{\link{det.eigen.Hermitian}}.
}
\examples{
x <- Matrix( sample(-3:3, size = 9, replace = T), nrow = 3, ncol = 3)
x[row(x) > col(x)] <- t(x)[row(x) > col(x)]  # construct symmetric matrix
class(x) <- Matrix.class(x)
det(x)
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
