\name{SVD}
\title{Singular Value Decomposition of a Matrix}
\usage{
SVD(x, nu = min(n,p), nv = min(n,p))
}
\alias{SVD}
\arguments{
  \item{x}{a matrix whose SVD decomposition is to be computed.}
  \item{nu}{the number of left eigenvectors to be computed.
    This must be one of \code{0}, \code{nrow(x)} and \code{ncol(x)}.}
  \item{nv}{the number of right eigenvectors to be computed.
    This must be one of \code{0}, and \code{ncol(x)}.}
}
\description{
  Compute the singular-value decomposition of a rectangular matrix.
}
\details{
  \code{SVD} provides an interface to the Lapack routine DGESVD.
  The singular value decomposition plays an important role in many
  statistical techniques.
}
\value{
  The SVD decomposition of the matrix as computed by Lapack,
  \deqn{ \bold{X = U D V'},} where \eqn{\bold{U}} and \eqn{\bold{V}} are
  orthogonal, \eqn{\bold{V'}} means \emph{V transposed}, and
  \eqn{\bold{D}} is a diagonal matrix with the singular
  values \eqn{D_{ii}}{D[i,i]}.  Equivalently, \eqn{\bold{D = U' X V}},
  which is verified in the examples, below.

  The components in the returned value correspond directly
  to the values returned by DGESVD.
  \item{d}{a vector containing the singular values of \code{x}.}
  \item{u}{a matrix whose columns contain the left eigenvectors of \code{x}.}
  \item{vt}{a matrix whose \bold{rows} contain the right eigenvectors of
      \code{x}.  That is, the \code{vt} component is the transpose of
      \eqn{\bold{V}} in the above decomposition.  This is different from
      the LINPACK result which returns \eqn{\bold{V}}.}
}
\references{
}
\seealso{\code{\link{svd}}, \code{\link{eigen}}, \code{\link{qr}}.}
\examples{
str(X <- hilbert(9)[,1:6])
str(s <- SVD(X))
Eps <- 100 * .Machine$double.eps

D <- diag(s$d)
stopifnot(abs(X - s$u \%*\% D \%*\% s$vt) < Eps)#  X = U D V'
stopifnot(abs(D - t(s$u) \%*\% X \%*\% t(s$vt)) < Eps)#  D = U' X V

X <- cbind(1, 1:7)
str(s <- SVD(X)); D <- diag(s$d)
stopifnot(abs(X - s$u \%*\% D \%*\% s$vt) < Eps)#  X = U D V'
stopifnot(abs(D - t(s$u) \%*\% X \%*\% t(s$vt)) < Eps)#  D = U' X V
}
\keyword{algebra}
\keyword{array}
