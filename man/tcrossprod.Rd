\name{tcrossprod}
\docType{methods}
\alias{tcrossprod-methods}
\alias{tcrossprod,cscMatrix-method}
\alias{tcrossprod}
\title{Cross-product of transpose}
\description{
  Take the cross-product of the transpose of a matrix.
  This is formally equivalent to, but faster than, the
  call \code{x \%*\% t(x)}.
}
\usage{
tcrossprod(x)
}
\arguments{
  \item{x}{a matrix-like object}
}
\details{
  For some classes in the \code{Matrix} package, such as the
  \code{\link{cscMatrix-class}}, it is much faster to calculate the
  cross-product of the transpose directly instead of calculating the
  transpose first and then its cross-product.
}
\value{
  An object of an appropriate symmetric matrix class.
}
\section{Methods}{
  \describe{
    \item{x = "cscMatrix"}{method for compressed, sparse,
      column-oriented matrices.}
  }
}
\seealso{\code{\link[base]{crossprod}}}
%\examples{}
\keyword{methods}
\keyword{algebra}
