\name{lu.Hermitian}
\alias{lu.Hermitian}
\title{
Triangular Decomposition of a Hermitian Matrix
}
\description{
Computes a symmetric indefinite factorization of a real symmetric
or complex Hermitian Matrix.
}
\usage{
lu.Hermitian(x, lower = T, norm.comp = T, tune = <see below>, 
             workspace = <see below>)
}
\arguments{
\item{x}{
a real symmetric or complex Hermitian matrix. Missing values (NAs) are not 
allowed.
}
\item{lower}{
a logical variable telling whether to use the lower or upper triangle of the 
matrix for the factorization. The default is to select the lower triangle.
}
\item{norm.comp}{
a logical variable telling whether or not to compute the one or infinity norm
of the matrix and return it as an attribute. The norm should be computed if
\code{solve} is going to eventually be applied to the factorization with condition
estimation. The default is to compute the norm.
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
A symmetric indefinite factorization of \code{x} is returned as
an object of class \code{c("lu.Hermitian", "decomp")}.
A copy of the call to \code{lu.Hermitian} is included as an 
attributes, as is the one or infinity norm of the matrix if requested 
(used for condition estimation).
}
\details{
Based on the functions dsytrf and zhetrf from LAPACK (Anderson et al. 1994).
For any symmetric (Hermitian) matrix \code{X}, there is a row permutation \code{P},
a triangular matrix \code{T} whose diagonal elements equal 1, and a block diagonal
matrix \code{B} with either 1 by 1 or symmetric 2 by 2 blocks on the diagonal, 
such that \code{P X t(P) = T B t(T)}. 
The matrix \code{T} is either lower or upper triangular (depending on which 
triangle is selected by \code{lower}).
This product is called a \code{symmetric indefinite factorization} of 
\code{X}. If all of the blocks of \code{B} are 1 by 1 (i. e. \code{B} is diagonal), and 
if all of the diagonal elements of \code{B} are positive, the the matrix is
positive definite and the factorization is a variant of the Cholesky
decomposition with pivoting.
}
\references{
Anderson, E., et al. (1994).
\emph{LAPACK User's Guide,}
2nd edition, SIAM, Philadelphia.
}
\seealso{
\code{\link{expand.lu.Hermitian}}, \code{\link{facmul.lu.Hermitian}}, \code{\link{.laenv}}, 
\code{\link{lu.Hermitian.object}},  \code{\link{lu.Matrix}}.
}
\examples{
 n <- 5
 x <- Matrix( rnorm(n*n), nrow = n, ncol = n)
 x[row(x) > col(x)] <- t(x)[row(x) > col(x)]  # construct symmetric matrix
 class(x) <- Matrix.class(x)
 lu(x)                           # symmetric indefinite factorization of A
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
