\name{facmul-methods}
\title{Multiplication by Factors from Matrix Factorizations}
%
\docType{methods}
\keyword{arith}
\keyword{array}
\keyword{methods}
%
\alias{facmul}
\alias{facmul-methods}
%
\description{
  Multiplies a matrix or vector on the left or right by a factor
  from a matrix factorization or its transpose.
}
\usage{
facmul(x, factor, y, trans = FALSE, left = TRUE, \dots)
}
\arguments{
  \item{x}{a \code{\linkS4class{MatrixFactorization}} object.}
  \item{factor}{a character string indicating a factor in the
    factorization represented by \code{x}, typically an element
    of \code{names(\link{expand2}(x, \dots))}.}
  \item{y}{a matrix or vector to be multiplied on the left or right
    by the factor or its transpose.}
  \item{trans}{a logical indicating if the transpose of the
    factor should be used, rather than the factor itself.}
  \item{left}{a logical indicating if the \code{y} should be
    multiplied on the left by the factor, rather than on the right.}
  \item{\dots}{further arguments passed to or from methods.}
}
\value{
  The value of \code{op(M) \%*\% y} or \code{y \%*\% op(M)},
  depending on \code{left}, where \code{M} is the factor
  (always \emph{without} \code{dimnames}) and \code{op(M)}
  is \code{M} or \code{t(M)}, depending on \code{trans}.
}
\details{
  \code{facmul} is experimental and currently no methods are
  exported from \pkg{Matrix}.
}
\examples{
\dontshow{ % for R_DEFAULT_PACKAGES=NULL
library(stats, pos = "package:base", verbose = FALSE)
}
## Conceptually, methods for 'facmul' _would_ behave as follows ...
\dontrun{
n <- 3L
x <- lu(Matrix(rnorm(n * n), n, n))
y <- rnorm(n)
L <- unname(expand2(x)[[nm <- "L"]])
stopifnot(exprs = {
    all.equal(facmul(x, nm, y, trans = FALSE, left =  TRUE), L \%*\% y)
    all.equal(facmul(x, nm, y, trans = FALSE, left = FALSE), y \%*\% L)
    all.equal(facmul(x, nm, y, trans =  TRUE, left =  TRUE),  crossprod(L, y))
    all.equal(facmul(x, nm, y, trans =  TRUE, left = FALSE), tcrossprod(y, L))
})
}
}
