\name{facmul.qr.Matrix}
\alias{facmul.qr.Matrix}
\title{
Multiplication by Factors from a QR Decomposition
}
\description{
Performs multiplication by any of the factors in a QR 
decomposition.
}
\usage{
facmul.qr.Matrix(x, factor, y, transpose = F, left = T, 
                   ncols, full = T, tune)
}
\arguments{
\item{x}{
an object of class \code{"qr.Matrix"}.
}
\item{factor}{
either "Q", "R", or "P",  where "Q" selects the orthonormal 
(unitary if complex) factor, "R" selects the triangular factor,
and "P selects the column permutation applied to the original matrix
in order to obtain the decomposition.
}
\item{y}{
a numeric or complex vector or matrix. The default is the identity 
matrix of the order equal to the number of rows (\code{left = T}) or 
number of columns (\code{left = F}) in the matrix underlying \code{x}, so that either 
the factor itself (\code{transpose = F}) or its transpose 
(\code{transpose = T}) are returned when \code{y} is missing.
}
\item{transpose}{
a logical value telling whether or not to multiply by the (if complex, 
conjugate) transpose of the factor. The default is multiplication by the 
factor itself.
}
\item{left}{
a logical value telling whether the multiplication should occur on the left 
of \code{y} or on the right. The default is multiplication on the left.
}
\item{ncols}{
gives the number of leading columns of the matrix underlying \code{x} 
(after pivoting) that are to be included in the orthogonal factor.
By default, the largest possible value is used for \code{ncols}, which is equal
to the minimum dimension of the matrix underlying \code{x}.
}
\item{full}{
if the matrix underlying \code{x} has row dimension m and column dimension n,
with m greater than n, then if \code{full = T}, Q is m x m and R is m x n, while 
if \code{full = F}, Q is m by n and R is n x n.
}
\item{tune}{
a integer vector or list of named tuning parameters that may affect
computational efficiency when no pivoting is done. The relevant parameters
are the blocksize parameters \code{NB}, \code{NBMIN}, and \code{NX} as described in \code{.laenv}.
The default is to use the settings in \code{.laenv}, which in S-Plus is initialized
to use an unblocked algorithm (\code{NB = 1}), but may be changed by the user.
}}
\value{
An object of class "Matrix" corresponding to the desired matrix product.
}
\details{
Based on the functions dormqr, dorgqr, zungqr, and zunmqr from LAPACK 
(Anderson et al. 1994).


The orthonormal/unitary factor in an \code{"qr.Matrix"} object is stored in a 
compact form representing elementary transformations whose product give the 
full factor. This function enables products involving this to be obtained 
efficiently without its explicit formation.
}
\references{
Anderson, E., et al. (1994).
\emph{LAPACK User's Guide,}
2nd edition, SIAM, Philadelphia.
}
\seealso{
\code{\link{qr.Matrix}}, \code{\link{expand.qr.Matrix}}
}
\examples{
 x <- Matrix( rnorm(36), nrow = 9, ncol = 4)
 z <- qr(x, pivot = T)                                # QR factorization of x
 max(abs(facmul(z,"Q", facmul(z,"R")) - x[,z$pivot])) # test product of Q and R
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
