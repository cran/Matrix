\name{eigen.Hermitian}
\alias{eigen.Hermitian}
\title{
Eigenvalue Decomposition of a Hermitian Matrix
}
\description{
Computes the eigenvalue decomposition of a real symmetric or complex
Hermitian Matrix.
}
\usage{
eigen.Hermitian(x, vectors=T, lower=T, value.range=<see below>, 
                index.range=<see below>, tol=<see below>,
                tune=<see below>, workspace=<see below>)
}
\arguments{
\item{x}{
a real symmetric or complex Hermitian matrix. 
Missing values (NAs) are not allowed.
}
\item{vectors}{
a logical value telling whether or not to compute eigenvectors.
The default is to compute all eigenvectors.
}
\item{lower}{
a logical variable telling whether the lower or upper triangle of \code{x} should 
be used for the computations. The default is to select the lower triangle.
}
\item{value.range}{
the range for eigenvalues. All eigenvalues on the half-open interval
\code{(min(value.range), max(value.range)]} will be found.
The default is to find all eigenvalues.
}
\item{index.range}{
the indeces of the eigenvalues of smallest and largest magnitude to be 
returned. The default is to find all eigenvalues.
}
\item{tol}{
absolute error tolerance for the eigenvalues. An approximate eigenvalue is
accepted as converged when it is determined to lie in an interval \code{[a,b]} of
width less than or equal to \code{tol + .Machine$double.eps * max(abs(a),abs(b))}.
If \code{tol} is non-positive, then the value \code{tol * .Machine$double.eps * T1},
where \code{T1} is the one norm of the tridiagonal reduction of \code{x}, is used.
The default is to use the latter definition of \code{tol}.
}
\item{tune}{
a integer vector or list of named tuning parameters that may affect
computational efficiency when no pivoting is done. The relevant parameters
are the blocksize parameters \code{NB}, \code{NBMIN}, and \code{NX} as described under
\code{.laenv}.
}
\item{workspace}{
workspace provided to the underlying software.
The default is to use the optimum value relative to the tuning parameters.
The optimal workspace for the problem (with the given tuning parameters) is
included as part of the output attributes.
}}
\value{
If no eigenvectors are requested then the result is a vector of eigenvalues.
Otherwise the result is a list with two named components \code{"values"} and
\code{"vectors"}.
If the entire eigenvalue decomposition is requested,
then the result is an object of class \code{c("eigen.Hermitian","decomp")}.
Workspace information and a copy of the call to \code{eigen.Hermitian} are included
as attributes.
}
\details{
Based on the functions dsyevx and zheevx from LAPACK (Anderson et al. 1994).
At most one of the parameters \code{value.range} and \code{index.range} may be specified
in any call.
Parameters in \code{tune} may affect performance and workspace requirements,
but the optimal settings are both machine and problem dependent.
Users may want to experiment with \code{tune} (see Anderson et al., p. 72-74) to
obtain improved performance.
}
\section{BUGS}{
The optimal workspace is sometimes computed incorrectly. 
A bug report has been sent to the Lapack maintainers.
}
\references{
Anderson, E., et al. (1994).
\emph{LAPACK User's Guide,}
2nd edition, SIAM, Philadelphia.
}
\seealso{
\code{\link{eigen.Hermitian.object}}, \code{\link{eigen.Matrix}}, \code{\link{.laenv}}
}
\examples{
n <- 5
x <- Matrix( rnorm(n*n), nrow = n, ncol = n)
x[row(x) < col(x)] <- t(x)[row(x) < col(x)]  # x is a real symmetric matrix
class(x) <-  Matrix.class(x)
eigen(x)                                     # eigenvalue decomposition of x 
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
