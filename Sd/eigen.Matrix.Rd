\name{eigen.Matrix}
\alias{eigen.Matrix}
\title{
Eigenvalue Decomposition of a Matrix
}
\description{
Computes the eigenvalue decomposition of a Matrix.
}
\usage{
eigen.Matrix(x, vectors = T, balance = F, rcond = F, schur = F, 
             tune = <see below>, workspace = <see below>)
}
\arguments{
\item{x}{
numeric or complex square Matrix inheriting from class \code{"Matrix"}.
Missing values (NAs) are not allowed.
}
\item{vectors}{
either a logical value telling whether or not to compute eigenvectors, or
else a logical vector of length 2 telling whether or not to compute the
left and right eigenvectors, respectively. The default is to compute both 
sets of eigenvectors.
}
\item{balance}{
either a logical value or a logical vector of length 2 telling whether or not
to perform permutations on the matrix to make it more nearly upper triangular 
and to perform diagonal scaling to make the rows and columns more equal in 
norm, respectively, prior to computing the eigenvalue decomposition.
Column permutation will not change the condition numbers 
(in exact arithmetic), but diagonal scaling will. 
Both may affect the performance of the algorithm.
By default, no balancing is performed.
}
\item{rcond}{
either a logical value or a logical vector of length 2 telling whether or 
not to compute reciprocal condition numbers for the eigenvalues and right 
eigenvectors, respectively, and return them as attributes. 
By default, neither is computed.
}
\item{schur}{
tells whether or not to return a matrix whose upper portion contains the 
Schur form (of the balanced matrix if balancing is requested) as an
attribute. By default, the Schur form is not returned.
}
\item{tune}{
a integer vector or list of named tuning parameters that may affect
computational efficiency. The relevant parameters are
\code{NB}, \code{NBMIN}, \code{NX}, \code{NS}, and \code{MAXB}, as described in \code{.laenv}.
The default is to use the settings in \code{.laenv}, which are initialized in
S-Plus but may be changed by the user.
}
\item{workspace}{
workspace provided to the underlying software.
The default is to use the minimum allowed value.
The optimal workspace for the problem (for the given \code{tune} specifications) is
included as part of the output attributes.
}}
\value{
The result is an object of class \code{c("eigen.Matrix", "decomp")}.
When the matrix is real, and there are complex eigenvalues, the
corresponding eigenvectors will be represented by their real and imaginary
parts, respectively, rather than by complex conjugate pairs.
They can be transformed into matrices of eigenvectors via the function
\code{expand}.
Attributes include the one-norm of the balanced matrix, balancing information
(if balancing is requested via \code{balance}), condition numbers (if requested via
\code{rcond}), a matrix containing the Schur form (if requested via \code{schur}),
workspace information, and a copy of the call to \code{eigen.Matrix}.
}
\details{
Based on the functions dgeevx and zgeevx from Lapack (Anderson et al. 1994).
Parameters in \code{tune} may affect performance and workspace requirements,
but the optimal settings are both machine and problem dependent.
Users may want to experiment with \code{tune} (see Anderson et al., p. 72-74) to
obtain improved performance.
}
\references{
Anderson, E., et al. (1994).
\emph{LAPACK User's Guide,}
2nd edition, SIAM, Philadelphia.
}
\seealso{
\code{\link{eigen.Hermitian}}, \code{\link{eigen.Matrix.object}}, \code{\link{expand.eigen.Matrix}},
\code{\link{.laenv}}, \code{\link{schur.Matrix}}, \code{\link{svd.Matrix}}.
}
\examples{
n <- 9
x <- Matrix( rnorm(n*n), nrow = n, ncol = n) 
eigen(x, vectors = F)$values              # eigenvalues of x
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
