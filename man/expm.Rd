\name{expm}
\alias{expm}
\alias{expm,dgeMatrix-method}
\title{Matrix exponential}
\description{
  Compute the exponential of a matrix.
}
\usage{
expm(x)
}
\arguments{
  \item{x}{an \R object that inherits from the \code{Matrix} class}
}
\details{
  The exponential of a matrix is defined as the infinite Taylor
  series \code{expm(A) = I + A + A^2/2! + A^3/3! + ...} (although this is
  definitely not the way to compute it).  The method for the
  \code{dgeMatrix} class uses Ward's diagonal Pade' approximation with
  three step preconditioning.
}
\value{
  The matrix exponential of \code{x}.
}
\references{}%% << FIXME
\author{This is a simple translation of the
  implementation of the corresponding Octave function contributed to the
  Octave project by A. Scottedward Hodel
  \email{A.S.Hodel@Eng.Auburn.EDU}.
}
%\note{}
\seealso{\code{\link{Schur}}}
\examples{
(m1 <- Matrix(c(1,0,1,1), nc = 2))
(e1 <- expm(m1)) ; e <- exp(1)
stopifnot(all.equal(e1@x, c(e,0,e,e), tol = 1e-15))
(m2 <- Matrix(c(-49, -64, 24, 31), nc = 2))
(e2 <- expm(m2))
(m3 <- Matrix(cbind(0,rbind(6*diag(3),0)), nc = 4))
(e3 <- expm(m3))
}
\keyword{algebra}
\keyword{math}
