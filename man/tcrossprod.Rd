\name{tcrossprod}
\docType{methods}
\alias{tcrossprod-methods}
\alias{tcrossprod,dgCMatrix-method}
\alias{tcrossprod,dgTMatrix-method}
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
  \code{\link{dgCMatrix-class}}, it is much faster to calculate the
  cross-product of the transpose directly instead of calculating the
  transpose first and then its cross-product.
}
\value{
  An object of an appropriate symmetric matrix class.
}
\section{Methods}{
  \describe{
    \item{x = "dgCMatrix"}{method for compressed, sparse,
      column-oriented matrices.}
  }
}
\seealso{\code{\link[base]{crossprod}}}
\examples{
 ## A random sparce "incidence" matrix :
 m <- matrix(0, 400, 500)
 set.seed(12)
 m[runif(314, 0, length(m))] <- 1
 mm <- as(m, "dgCMatrix")
 object.size(m) / object.size(mm) # smaller by a factor of 242.88

 ## tcrossprod() is very fast:
 system.time(tCmm <- tcrossprod(mm))# "0" practically
 system.time(cm <- crossprod(t(m))) #

 stopifnot(identical(cm, as(tCmm, "matrix")))

 ## show "sparse" (sub) matrix
 tc <- cm[1:16, 1:30]
 storage.mode(tc) <- "character" ; tc[tc == "0"] <- "."
 dimnames(tc)[[2]] <- rep("", ncol(tc))
 noquote(tc)
}
\keyword{methods}
\keyword{algebra}
