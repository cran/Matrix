\name{BunchKaufman-methods}
\title{Methods for Bunch-Kaufman Factorization}
%
\docType{methods}
\keyword{algebra}
\keyword{array}
\keyword{methods}
%
\alias{BunchKaufman}
\alias{BunchKaufman-methods}
%
\alias{BunchKaufman,dspMatrix-method}
\alias{BunchKaufman,dsyMatrix-method}
\alias{BunchKaufman,matrix-method}
%
\description{
  Computes the Bunch-Kaufman factorization of an \eqn{n \times n}{n-by-n}
  real, symmetric matrix \eqn{A}, which has the general form
  \deqn{A = U D_{U} U' = L D_{L} L'}{A = U * DU * U' = L * DL * L'}
  where
  \eqn{D_{U}}{DU} and \eqn{D_{L}}{DL} are symmetric, block diagonal
  matrices composed of \eqn{b_{U}}{bU} and \eqn{b_{L}}{bL}
  \eqn{1 \times 1}{1-by-1} or \eqn{2 \times 2}{2-by-2} diagonal blocks;
  \eqn{U = \prod_{k = 1}^{b_{U}} P_{k} U_{k}}{U = prod(Pk * Uk : k = 1,...,bU)}
  is the product of \eqn{b_{U}}{bU} row-permuted unit upper triangular
  matrices, each having nonzero entries above the diagonal in 1 or 2 columns;
  and
  \eqn{L = \prod_{k = 1}^{b_{L}} P_{k} L_{k}}{L = prod(Pk * Lk : k = 1,...,bL)}
  is the product of \eqn{b_{L}}{bL} row-permuted unit lower triangular
  matrices, each having nonzero entries below the diagonal in 1 or 2 columns.
  
  Methods are built on LAPACK routines \code{dsytrf} and \code{dsptrf}.
}
\usage{
BunchKaufman(x, \dots)
\S4method{BunchKaufman}{dsyMatrix}(x, warnSing = TRUE, \dots)
\S4method{BunchKaufman}{dspMatrix}(x, warnSing = TRUE, \dots)
\S4method{BunchKaufman}{matrix}(x, uplo = "U", \dots)
}
\arguments{
  \item{x}{a \link[=is.finite]{finite} symmetric matrix or
    \code{\linkS4class{Matrix}} to be factorized.
    If \code{x} is square but not symmetric, then it will be
    \emph{treated} as symmetric; see \code{uplo}.}
  \item{warnSing}{a logical indicating if a \link{warning} should
    be signaled for singular \code{x}.}
  \item{uplo}{a string, either \code{"U"} or \code{"L"},
    indicating which triangle of \code{x} should be used
    to compute the factorization.}
  \item{\dots}{further arguments passed to or from methods.}
}
\value{
  An object representing the factorization, inheriting from
  virtual class \code{\linkS4class{BunchKaufmanFactorization}}.
  The specific class is \code{\linkS4class{BunchKaufman}} unless
  \code{x} inherits from virtual class \code{\linkS4class{packedMatrix}},
  in which case it is \code{\linkS4class{pBunchKaufman}}.
}
\seealso{
  Classes \code{\linkS4class{BunchKaufman}} and
  \code{\linkS4class{pBunchKaufman}} and their methods.
  
  Classes \code{\linkS4class{dsyMatrix}} and
  \code{\linkS4class{dspMatrix}}.
  
  Generic functions \code{\link{expand1}} and \code{\link{expand2}},
  for constructing matrix factors from the result.
  
  Generic functions \code{\link{Cholesky}}, \code{\link{Schur}},
  \code{\link{lu}}, and \code{\link{qr}},
  for computing other factorizations.
}
\references{
  The LAPACK source code, including documentation; see
  \url{https://netlib.org/lapack/double/dsytrf.f} and
  \url{https://netlib.org/lapack/double/dsptrf.f}.

  Golub, G. H., & Van Loan, C. F. (2013).
  \emph{Matrix computations} (4th ed.).
  Johns Hopkins University Press.
  \doi{10.56021/9781421407944}
}
\examples{
\dontshow{ % for R_DEFAULT_PACKAGES=NULL
library(utils, pos = "package:base", verbose = FALSE)
}
showMethods("BunchKaufman", inherited = FALSE)
set.seed(0)

data(CAex, package = "Matrix")
class(CAex) # dgCMatrix
isSymmetric(CAex) # symmetric, but not formally

A <- as(CAex, "symmetricMatrix")
class(A) # dsCMatrix

## Have methods for denseMatrix (unpacked and packed),
## but not yet sparseMatrix ...
\dontrun{
(bk.A <- BunchKaufman(A))
}
(bk.A <- BunchKaufman(as(A, "unpackedMatrix")))

## A ~ U DU U' in floating point
str(e.bk.A <- expand2(bk.A), max.level = 2L)
stopifnot(all.equal(as(A, "matrix"), as(Reduce(`\%*\%`, e.bk.A), "matrix")))
}
