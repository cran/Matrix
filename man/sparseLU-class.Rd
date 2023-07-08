\name{sparseLU-class}
\title{Sparse LU Factorizations}
%
\docType{class}
\keyword{algebra}
\keyword{array}
\keyword{classes}
%
\alias{sparseLU-class}
%
\alias{determinant,sparseLU,logical-method}
%
\description{
  \code{sparseLU} is the class of sparse, row- and column-pivoted
  LU factorizations of \eqn{n \times n}{n-by-n} real matrices \eqn{A},
  having the general form
  \deqn{P_{1} A P_{2} = L U}{P1 * A * P2 = L * U}
  or (equivalently)
  \deqn{A = P_{1}' L U P_{2}'}{A = P1' * L * U * P2'}
  where
  \eqn{P_{1}}{P1} and \eqn{P_{2}}{P2} are permutation matrices,
  \eqn{L} is a unit lower triangular matrix, and
  \eqn{U} is an upper triangular matrix.
}
\section{Slots}{
  \describe{
    \item{\code{Dim}, \code{Dimnames}}{inherited from virtual class
      \code{\linkS4class{MatrixFactorization}}.}
    \item{\code{L}}{an object of class \code{\linkS4class{dtCMatrix}},
      the unit lower triangular \eqn{L} factor.}
    \item{\code{U}}{an object of class \code{\linkS4class{dtCMatrix}},
      the upper triangular \eqn{U} factor.}
    \item{\code{p}, \code{q}}{0-based integer vectors of length
      \code{Dim[1]},
      specifying the permutations applied to the rows and columns of
      the factorized matrix.  \code{q} of length 0 is valid and
      equivalent to the identity permutation, implying no column pivoting.
      Using \R{} syntax, the matrix \eqn{P_{1} A P_{2}}{P1 * A * P2}
      is precisely \code{A[p+1, q+1]}
      (\code{A[p+1, ]} when \code{q} has length 0).}
  }
}
\section{Extends}{
  Class \code{\linkS4class{LU}}, directly.
  Class \code{\linkS4class{MatrixFactorization}}, by class
  \code{\linkS4class{LU}}, distance 2.
}
\section{Instantiation}{
  Objects can be generated directly by calls of the form
  \code{new("sparseLU", ...)}, but they are more typically obtained
  as the value of \code{\link{lu}(x)} for \code{x} inheriting from
  \code{\linkS4class{sparseMatrix}} (often \code{\linkS4class{dgCMatrix}}).
}
\section{Methods}{
  \describe{
    \item{\code{determinant}}{\code{signature(from = "sparseLU", logarithm = "logical")}:
      computes the determinant of the factorized matrix \eqn{A}
      or its logarithm.}
    \item{\code{expand}}{\code{signature(x = "sparseLU")}:
      see \code{\link{expand-methods}}.}
    \item{\code{expand1}}{\code{signature(x = "sparseLU")}:
      see \code{\link{expand1-methods}}.}
    \item{\code{expand2}}{\code{signature(x = "sparseLU")}:
      see \code{\link{expand2-methods}}.}
    \item{\code{solve}}{\code{signature(a = "sparseLU", b = .)}:
      see \code{\link{solve-methods}}.}
  }
}
\seealso{
  Class \code{\linkS4class{denseLU}} for dense LU factorizations.
  
  Class \code{\linkS4class{dgCMatrix}}.

  Generic functions \code{\link{lu}},
  \code{\link{expand1}} and \code{\link{expand2}}.
}
\references{
  Davis, T. A. (2006).
  \emph{Direct methods for sparse linear systems}.
  Society for Industrial and Applied Mathematics.
  \doi{10.1137/1.9780898718881}

  Golub, G. H., & Van Loan, C. F. (2013).
  \emph{Matrix computations} (4th ed.).
  Johns Hopkins University Press.
  \doi{10.56021/9781421407944}
}
\examples{
\dontshow{ % for R_DEFAULT_PACKAGES=NULL
library(stats, pos = "package:base", verbose = FALSE)
library(utils, pos = "package:base", verbose = FALSE)
}
showClass("sparseLU")
set.seed(2)

A <- as(readMM(system.file("external", "pores_1.mtx", package = "Matrix")),
        "CsparseMatrix")
(n <- A@Dim[1L])

## With dimnames, to see that they are propagated :
dimnames(A) <- dn <- list(paste0("r", seq_len(n)),
                          paste0("c", seq_len(n)))

(lu.A <- lu(A))
str(e.lu.A <- expand2(lu.A), max.level = 2L)

ae1 <- function(a, b, ...) all.equal(as(a, "matrix"), as(b, "matrix"), ...)
ae2 <- function(a, b, ...) ae1(unname(a), unname(b), ...)

## A ~ P1' L U P2' in floating point
stopifnot(exprs = {
    identical(names(e.lu.A), c("P1.", "L", "U", "P2."))
    identical(e.lu.A[["P1."]],
              new("pMatrix", Dim = c(n, n), Dimnames = c(dn[1L], list(NULL)),
                  margin = 1L, perm = invertPerm(lu.A@p, 0L, 1L)))
    identical(e.lu.A[["P2."]],
              new("pMatrix", Dim = c(n, n), Dimnames = c(list(NULL), dn[2L]),
                  margin = 2L, perm = invertPerm(lu.A@q, 0L, 1L)))
    identical(e.lu.A[["L"]], lu.A@L)
    identical(e.lu.A[["U"]], lu.A@U)
    ae1(A, with(e.lu.A, P1. \%*\% L \%*\% U \%*\% P2.))
    ae2(A[lu.A@p + 1L, lu.A@q + 1L], with(e.lu.A, L \%*\% U))
})

## Factorization handled as factorized matrix
b <- rnorm(n)
stopifnot(identical(det(A), det(lu.A)),
          identical(solve(A, b), solve(lu.A, b)))
}
