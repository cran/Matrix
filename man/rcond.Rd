\name{rcond}
\title{Estimate the Reciprocal Condition Number}
\alias{rcond}
% most methods are documented in <foo>Matrix-class.Rd
\alias{rcond,ANY,missing-method}
\alias{rcond,matrix,character-method}
\alias{rcond,Matrix,character-method}
\alias{rcond,ldenseMatrix,character-method}
\alias{rcond,ndenseMatrix,character-method}
%
\usage{
rcond(x, type, \dots)
}
\description{
  Estimate the reciprocal of the condition number of a matrix.

  This is a generic function with several methods, as seen by
  \code{\link{showMethods}(rcond)}.
}
\arguments{
  \item{x}{an \R object that inherits from the \code{Matrix} class.}
  \item{type}{
    Character indicating the type of norm to be used in the estimate.
    The default is \code{"O"} for the 1-norm.  The other possible value is
    \code{"I"} for the infinity norm, see also \code{\link{norm}}.
  }
  \item{\dots}{further arguments passed to or from other methods.}
}
\value{
  An estimate of the reciprocal condition number of \code{x}.
}
\section{BACKGROUND}{
  The condition number of a matrix is the product of the
  \code{\link{norm}} of the matrix and the norm of its inverse (or
  pseudo-inverse).
  The condition number takes on values between 1 and infinity,
  inclusive, and can be viewed as a factor by which errors in solving
  linear systems with this matrix as coefficient matrix could be
  magnified.

  \code{rcond()} computes the \emph{reciprocal} condition number with
  values in \eqn{[0,1]} and can be viewed as a scaled measure of how
  close a matrix is to being rank deficient (aka \dQuote{singular}).

  Condition numbers are usually estimated, since exact computation is
  costly in terms of floating-point operations.  An (over) estimate of
  reciprocal condition number is given, since by doing so overflow is
  avoided.  Matrices are well-conditioned if the reciprocal condition
  number is near 1 and ill-conditioned if it is near zero.
}
\seealso{
  \code{\link{norm}}, \code{\link[base]{solve}}.
}
\references{
  Golub, G., and Van Loan, C. F. (1989).
  \emph{Matrix Computations,}
  2nd edition, Johns Hopkins, Baltimore.
}
\examples{
x <- Matrix(rnorm(9), 3, 3)
rcond(x)
## typically "the same" (with more computational effort):
1 / (norm(x) * norm(solve(x)))

rcond(Hilbert(9))  # should be about 9.1e-13
}
\keyword{array}
\keyword{algebra}
