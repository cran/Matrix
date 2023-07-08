\name{chol-methods}
\title{Compute the Cholesky Factor of a Matrix}
\docType{methods}
%
\keyword{algebra}
\keyword{array}
\keyword{methods}
%
\alias{chol}
\alias{chol-methods}
%
\alias{chol,ddiMatrix-method}
\alias{chol,diagonalMatrix-method}
\alias{chol,dsCMatrix-method}
\alias{chol,dsRMatrix-method}
\alias{chol,dsTMatrix-method}
\alias{chol,dspMatrix-method}
\alias{chol,dsyMatrix-method}
\alias{chol,generalMatrix-method}
\alias{chol,symmetricMatrix-method}
\alias{chol,triangularMatrix-method}
%
\description{
  Computes the upper triangular Cholesky factor of an
  \eqn{n \times n}{n-by-n} real, symmetric, positive semidefinite
  matrix \eqn{A}, optionally after pivoting.
  That is the factor \eqn{L'} in
  \deqn{P_{1} A P_{1}' = L L'}{P1 * A * P1' = L * L'}
  or (equivalently)
  \deqn{A = P_{1}' L L' P_{1}}{A = P1' * L * L' * P1}
  where
  \eqn{P_{1}}{P1} is a permutation matrix.

  Methods for \code{\linkS4class{denseMatrix}} are built on
  LAPACK routines \code{dpstrf}, \code{dpotrf}, and \code{dpptrf},
  The latter two do not permute rows or columns,
  so that \eqn{P_{1}}{P1} is an identity matrix.
  
  Methods for \code{\linkS4class{sparseMatrix}} are built on
  CHOLMOD routines \code{cholmod_analyze} and \code{cholmod_factorize_p}.
}
\usage{
chol(x, \dots)
\S4method{chol}{dsyMatrix}(x, pivot = FALSE, tol = -1, \dots)
\S4method{chol}{dspMatrix}(x, \dots)
\S4method{chol}{dsCMatrix}(x, pivot = FALSE, \dots)
\S4method{chol}{ddiMatrix}(x, \dots)
\S4method{chol}{generalMatrix}(x, uplo = "U", \dots)
\S4method{chol}{triangularMatrix}(x, uplo = "U", \dots)
}
\arguments{
  \item{x}{a \link[=is.finite]{finite}, symmetric, positive
    semidefinite matrix or \code{\linkS4class{Matrix}} to
    be factorized.  If \code{x} is square but not symmetric,
    then it will be \emph{treated} as symmetric; see \code{uplo}.
    Methods for dense \code{x} require positive definiteness
    when \code{pivot = FALSE}.
    Methods for sparse (but not diagonal) \code{x} require
    positive definiteness unconditionally.}
  \item{pivot}{a logical indicating if the rows and columns
    of \eqn{x} should be pivoted.  Methods for sparse \code{x}
    employ the approximate minimum degree (AMD) algorithm
    in order to reduce fill-in, i.e., without regard for
    numerical stability.}
  \item{tol}{a \link[=is.finite]{finite} numeric tolerance,
    used only if \code{pivot = TRUE}.
    The factorization algorithm stops if the pivot is less than
    or equal to \code{tol}.  Negative \code{tol} is equivalent
    to \code{nrow(x) * .Machine$double.eps * max(diag(x))}.}
  \item{uplo}{a string, either \code{"U"} or \code{"L"},
    indicating which triangle of \code{x} should be used
    to compute the factorization.  The default is \code{"U"},
    even for lower triangular \code{x}, to be consistent with
    \code{\link[base]{chol}} from \pkg{base}.}
  \item{\dots}{further arguments passed to or from methods.}
}
\value{
  A matrix, \code{\linkS4class{triangularMatrix}},
  or \code{\linkS4class{diagonalMatrix}} representing
  the upper triangular Cholesky factor \eqn{L'}.
  The result is a traditional matrix if \code{x} is a
  traditional matrix, dense if \code{x} is dense, and
  sparse if \code{x} is sparse.
}
\details{
  For \code{x} inheriting from \code{\linkS4class{diagonalMatrix}},
  the diagonal result is computed directly and without pivoting,
  i.e., bypassing CHOLMOD.
  
  For all other \code{x}, \code{chol(x, pivot = value)} calls
  \code{\link{Cholesky}(x, perm = value, \dots)} under the hood.
  If you must know the permutation \eqn{P_{1}}{P1} in addition
  to the Cholesky factor \eqn{L'}, then call \code{\link{Cholesky}}
  directly, as the result of \code{chol(x, pivot = TRUE)} specifies
  \eqn{L'} but not \eqn{P_{1}}{P1}.
}
\seealso{
  The default method from \pkg{base}, \code{\link[base]{chol}},
  called for traditional matrices \code{x}.

  Generic function \code{\link{Cholesky}}, for more flexibility
  notably when computing the Cholesky \emph{factorization} and
  not only the \emph{factor} \eqn{L'}.
}
\references{
  The LAPACK source code, including documentation; see
  \url{https://netlib.org/lapack/double/dpstrf.f},
  \url{https://netlib.org/lapack/double/dpotrf.f}, and
  \url{https://netlib.org/lapack/double/dpptrf.f}.

  The CHOLMOD source code; see
  \url{https://github.com/DrTimothyAldenDavis/SuiteSparse},
  notably the header file \file{CHOLMOD/Include/cholmod.h}
  defining \code{cholmod_factor_struct}.

  Chen, Y., Davis, T. A., Hager, W. W., & Rajamanickam, S. (2008).
  Algorithm 887: CHOLMOD, supernodal sparse Cholesky factorization
  and update/downdate.
  \emph{ACM Transactions on Mathematical Software},
  \emph{35}(3), Article 22, 1-14.
  \doi{10.1145/1391989.1391995}

  Amestoy, P. R., Davis, T. A., & Duff, I. S. (2004).
  Algorithm 837: AMD, an approximate minimum degree ordering algorithm.
  \emph{ACM Transactions on Mathematical Software},
  \emph{17}(4), 886-905.
  \doi{10.1145/1024074.1024081}

  Golub, G. H., & Van Loan, C. F. (2013).
  \emph{Matrix computations} (4th ed.).
  Johns Hopkins University Press.
  \doi{10.56021/9781421407944}
}
\examples{
\dontshow{ % for R_DEFAULT_PACKAGES=NULL
library(stats, pos = "package:base", verbose = FALSE)
}
showMethods("chol", inherited = FALSE)
set.seed(0)

## ---- Dense ----------------------------------------------------------

## chol(x, pivot = value) wrapping Cholesky(x, perm = value)
selectMethod("chol", "dsyMatrix")

## Except in packed cases where pivoting is not yet available
selectMethod("chol", "dspMatrix")

## .... Positive definite ..............................................

(A1 <- new("dsyMatrix", Dim = c(2L, 2L), x = c(1, 2, 2, 5)))
(R1.nopivot <- chol(A1))
(R1 <- chol(A1, pivot = TRUE))

## In 2-by-2 cases, we know that the permutation is 1:2 or 2:1,
## even if in general 'chol' does not say ...

stopifnot(exprs = {
   all.equal(  A1           , as(crossprod(R1.nopivot), "dsyMatrix"))
   all.equal(t(A1[2:1, 2:1]), as(crossprod(R1        ), "dsyMatrix"))
   identical(Cholesky(A1)@perm, 2:1) # because 5 > 1
})

## .... Positive semidefinite but not positive definite ................

(A2 <- new("dpoMatrix", Dim = c(2L, 2L), x = c(1, 2, 2, 4)))
try(R2.nopivot <- chol(A2)) # fails as not positive definite
(R2 <- chol(A2, pivot = TRUE)) # returns, with a warning and ...

stopifnot(exprs = {
   all.equal(t(A2[2:1, 2:1]), as(crossprod(R2), "dsyMatrix"))
   identical(Cholesky(A2)@perm, 2:1) # because 4 > 1
})

## .... Not positive semidefinite ......................................

(A3 <- new("dsyMatrix", Dim = c(2L, 2L), x = c(1, 2, 2, 3)))
try(R3.nopivot <- chol(A3)) # fails as not positive definite
(R3 <- chol(A3, pivot = TRUE)) # returns, with a warning and ...

## _Not_ equal: see details and examples in help("Cholesky")
all.equal(t(A3[2:1, 2:1]), as(crossprod(R3), "dsyMatrix"))

## ---- Sparse ---------------------------------------------------------

## chol(x, pivot = value) wrapping
## Cholesky(x, perm = value, LDL = FALSE, super = FALSE)
selectMethod("chol", "dsCMatrix")

## Except in diagonal cases which are handled "directly"
selectMethod("chol", "ddiMatrix")

(A4 <- toeplitz(as(c(10, 0, 1, 0, 3), "sparseVector")))
(ch.A4.nopivot <- Cholesky(A4, perm = FALSE, LDL = FALSE, super = FALSE))
(ch.A4 <- Cholesky(A4, perm = TRUE, LDL = FALSE, super = FALSE))
(R4.nopivot <- chol(A4))
(R4 <- chol(A4, pivot = TRUE))

det4 <- det(A4)
b4 <- rnorm(5L)
x4 <- solve(A4, b4)

stopifnot(exprs = {
    identical(R4.nopivot, expand1(ch.A4.nopivot, "L."))
    identical(R4, expand1(ch.A4, "L."))
    all.equal(A4, crossprod(R4.nopivot))
    all.equal(A4[ch.A4@perm + 1L, ch.A4@perm + 1L], crossprod(R4))
    all.equal(diag(R4.nopivot), sqrt(diag(ch.A4.nopivot)))
    all.equal(diag(R4), sqrt(diag(ch.A4)))
    all.equal(sqrt(det4), det(R4.nopivot))
    all.equal(sqrt(det4), det(R4))
    all.equal(det4, det(ch.A4.nopivot, sqrt = FALSE))
    all.equal(det4, det(ch.A4, sqrt = FALSE))
    all.equal(x4, solve(R4.nopivot, solve(t(R4.nopivot), b4)))
    all.equal(x4, solve(ch.A4.nopivot, b4))
    all.equal(x4, solve(ch.A4, b4))
})
}
