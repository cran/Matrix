\name{solve.Matrix}
\title{Solve a System of Equations}
\usage{
\method{solve}{Matrix}(a, b, tol = 0, transpose = FALSE, \dots)
}
\alias{solve.Matrix}
\arguments{
  \item{a}{a numeric matrix containing the coefficients of the linear
    system.}
  \item{b}{a numeric vector or matrix giving the right-hand side(s) of
    the linear system.  If omitted, \code{b} is taken to be an identity
    matrix and \code{solve} will return the inverse of \code{a}.}
  \item{tol}{the tolerance for the reciprocal condition estimate.  If
      \code{tol} is negative, no condition estimation is done.
      Otherwise the reciprocal one norm condition estimate is computed 
      and the solution is attempted only when it is greater than \code{tol}.}
  \item{transpose}{a logical value.  If \code{TRUE}, the system to be
      solved is defined by the transpose of \code{a}.  The default is
      \code{FALSE}, for which the system to be solved is defined by
       \code{a}.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\description{
    This function solves the equation \code{a \%*\% x = b} or
    \code{t(a) \%*\% x = b} for \code{x}, where \code{b} can be either a
    vector or a matrix. 
}
\seealso{
  \code{\link{rcond.Matrix}}
}
\examples{
h8 <- hilbert(8); round(h8, 4)
solve(h8)
solve(h8, diag(8)) # similar result but not necessarily Hermitian
Hermitian.test(solve(h8))
Hermitian.test(solve(h8, diag(8)))
}
\keyword{algebra}
