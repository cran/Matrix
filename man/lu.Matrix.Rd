\name{lu.Matrix}
\alias{lu.Matrix}
\title{
  Triangular (LU) Decomposition of a Matrix
}
\description{
  Computes the LU decomposition of a matrix.
}
\usage{
lu.Matrix(x, norm.comp = T)
}
\arguments{
  \item{x}{
    numeric or complex matrix (does not have to be square). Missing values (NAs) 
    are not allowed.
  }
  \item{norm.comp}{
    a logical vector of length two telling whether or not to compute the 
    one and/or infinity norm of the matrix, respectively. These norms
    should be computed if \code{solve} is going to eventually be applied
    to the factorization with condition estimation. The infinity norm is
    needed for solves involving the underlying matrix, and the one norm
    is needed for solves involving its transpose. 
  }
}
\value{
  The LU decomposition of \code{x} is returned as an object of class
  \code{"lu.Matrix"}, which is a list with components \code{"l"}, a unit
  lower triangular matrix, \code{"u"}, an upper triangular matrix, and
  \code{"permutation"}, a row permutation matrix.  The one or infinity
  norm is included as an attribute, if requested.
  and/or infinity norms of the matrix if requested (used for condition
  estimation). 
}
\details{
  This function uses the Lapack subroutine DGETRF to form an LU
  factorization of a matrix X.
  The factorization has the form
  X = P * L * U
  where  P  is  a  permutation matrix, L is lower triangular
  with unit diagonal elements (lower trapezoidal if m >  n),
  and U is upper triangular (upper trapezoidal if m < n).
}
\references{
  Anderson, E., et al. (1994).
  \emph{LAPACK User's Guide,}
  2nd edition, SIAM, Philadelphia.
}
\examples{
x <- Matrix((1:16)^2, nrow = 4, ncol = 4)
lu.Matrix(x)
}
\keyword{algebra}
