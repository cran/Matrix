\name{lu.Matrix}
\alias{lu.Matrix}
\title{
Triangular (LU) Decomposition of a Matrix
}
\description{
Computes the LU decomposition of a matrix.
}
\usage{
lu.Matrix(x, norm.comp = T, tune = <see below>)
}
\arguments{
\item{x}{
numeric or complex matrix (does not have to be square). Missing values (NAs) 
are not allowed.
}
\item{norm.comp}{
a logical vector of length two telling whether or not to compute the one and/or
infinity norm of the matrix, respectively. These norms should be computed if
\code{solve} is going to eventually be applied to the factorization with condition
estimation. The infinity norm is needed for solves involving the underlying
matrix, and the one norm is needed for solves involving its transpose.
}
\item{tune}{
a integer vector or list of named tuning parameters that may affect
computational efficiency. The relevant parameter is the blocksize parameter
\code{NB} as described in \code{.laenv}.
}}
\value{
The LU decomposition of \code{x} is returned as
an object of class \code{c("lu.Matrix", "decomp")}.
A copy of the call to \code{lu.Matrix} is included as an attribute,
as well as the one and/or infinity norms of the matrix if requested
(used for condition estimation).
}
\details{
Based on the functions dgetrf and zgetrf from LAPACK (Anderson et al. 1994).
For any square matrix \code{X}, there is a row permutation \code{P}, 
a lower triangular matrix \code{L} whose diagonal elements equal 1, and an
upper triangular matrix \code{U}, such that \code{P X} = \code{L U}.
If \code{X} is not square, a similar factorization exists, with either \code{L} or \code{U} 
trapezoidal rather than triangular.
This product is called an \code{LU} decomposition of \code{X}.
The blocksize parameter (NB) in \code{tune} may affect performance, 
but the optimal settings are both machine and problem dependent.
Users may want to experiment with \code{tune} (see Anderson et al., p. 72-74)
to obtain improved performance.
}
\references{
Anderson, E., et al. (1994).
\emph{LAPACK User's Guide,}
2nd edition, SIAM, Philadelphia.
}
\seealso{
\code{\link{expand.lu.Matrix}}, \code{\link{facmul.lu.Matrix}}, \code{\link{.laenv}}, \code{\link{lu.Matrix.object}}, 
\code{\link{lu.Hermitian}}.
}
\examples{
x <- Matrix(rnorm(25), nrow = 5, ncol = 5)
lu.Matrix(x)
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
