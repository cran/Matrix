\name{denseLU-class}
\title{Dense LU Factorizations}
%
\docType{class}
\keyword{algebra}
\keyword{array}
\keyword{classes}
%
\alias{denseLU-class}
%
\alias{coerce,denseLU,dgeMatrix-method}
\alias{determinant,denseLU,logical-method}
%
\description{
  \code{denseLU} is the class of dense, row-pivoted LU factorizations
  of \eqn{m \times n}{m-by-n} real matrices \eqn{A},
  having the general form
  \deqn{P_{1} A = L U}{P1 * A = L * U}
  or (equivalently)
  \deqn{A = P_{1}' L U}{A = P1' * L * U}
  where
  \eqn{P_{1}}{P1} is an \eqn{m \times m}{m-by-m} permutation matrix, 
  \eqn{L} is an \eqn{m \times \text{min}(m,n)}{m-by-min(m,n)}
  unit lower trapezoidal matrix, and
  \eqn{U} is a  \eqn{\text{min}(m,n) \times n}{min(m,n)-by-n}
  upper trapezoidal matrix.  If \eqn{m = n}, then the factors
  \eqn{L} and \eqn{U} are triangular.
}
\section{Slots}{
  \describe{
    \item{\code{Dim}, \code{Dimnames}}{inherited from virtual class
      \code{\linkS4class{MatrixFactorization}}.}
    \item{\code{x}}{a numeric vector of length \code{prod(Dim)} storing
      the triangular \eqn{L} and \eqn{U} factors together in a packed
      format.  The details of the representation are specified by the
      manual for LAPACK routine \code{dgetrf}.}
    \item{\code{perm}}{an integer vector of length \code{min(Dim)}
      specifying the permutation \eqn{P_{1}}{P1} as a product of
      transpositions.   The corresponding permutation vector can
      be obtained as \code{\link{asPerm}(perm)}.}
  }
}
\section{Extends}{
  Class \code{\linkS4class{LU}}, directly.
  Class \code{\linkS4class{MatrixFactorization}}, by class
  \code{\linkS4class{LU}}, distance 2.
}
\section{Instantiation}{
  Objects can be generated directly by calls of the form
  \code{new("denseLU", ...)}, but they are more typically obtained
  as the value of \code{\link{lu}(x)} for \code{x} inheriting from
  \code{\linkS4class{denseMatrix}} (often \code{\linkS4class{dgeMatrix}}).
}
\section{Methods}{
  \describe{
    \item{\code{coerce}}{\code{signature(from = "denseLU", to = "dgeMatrix")}:
      returns a \code{\linkS4class{dgeMatrix}} with the dimensions
      of the factorized matrix \eqn{A}, equal to \eqn{L} below the
      diagonal and equal to \eqn{U} on and above the diagonal.}
    \item{\code{determinant}}{\code{signature(from = "denseLU", logarithm = "logical")}:
      computes the determinant of the factorized matrix \eqn{A}
      or its logarithm.}
    \item{\code{expand}}{\code{signature(x = "denseLU")}:
      see \code{\link{expand-methods}}.}
    \item{\code{expand1}}{\code{signature(x = "denseLU")}:
      see \code{\link{expand1-methods}}.}
    \item{\code{expand2}}{\code{signature(x = "denseLU")}:
      see \code{\link{expand2-methods}}.}
    \item{\code{solve}}{\code{signature(a = "denseLU", b = "missing")}:
      see \code{\link{solve-methods}}.}
  }
}
\seealso{
  Class \code{\linkS4class{sparseLU}} for sparse LU factorizations.
  
  Class \code{\linkS4class{dgeMatrix}}.

  Generic functions \code{\link{lu}},
  \code{\link{expand1}} and \code{\link{expand2}}.
}
\references{
  The LAPACK source code, including documentation; see
  \url{https://netlib.org/lapack/double/dgetrf.f}.
  
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
showClass("denseLU")
set.seed(1)

n <- 3L
(A <- Matrix(round(rnorm(n * n), 2L), n, n))

## With dimnames, to see that they are propagated :
dimnames(A) <- dn <- list(paste0("r", seq_len(n)),
                          paste0("c", seq_len(n)))

(lu.A <- lu(A))
str(e.lu.A <- expand2(lu.A), max.level = 2L)

## Underlying LAPACK representation
(m.lu.A <- as(lu.A, "dgeMatrix")) # which is L and U interlaced
stopifnot(identical(as(m.lu.A, "matrix"), `dim<-`(lu.A@x, lu.A@Dim)))

ae1 <- function(a, b, ...) all.equal(as(a, "matrix"), as(b, "matrix"), ...)
ae2 <- function(a, b, ...) ae1(unname(a), unname(b), ...)

## A ~ P1' L U in floating point
stopifnot(exprs = {
    identical(names(e.lu.A), c("P1.", "L", "U"))
    identical(e.lu.A[["P1."]],
              new(  "pMatrix", Dim = c(n, n), Dimnames = c(dn[1L], list(NULL)),
                  margin = 1L, perm = invertPerm(asPerm(lu.A@perm))))
    identical(e.lu.A[["L"]],
              new("dtrMatrix", Dim = c(n, n), Dimnames = list(NULL, NULL),
                  uplo = "L", diag = "U", x = lu.A@x))
    identical(e.lu.A[["U"]],
              new("dtrMatrix", Dim = c(n, n), Dimnames = c(list(NULL), dn[2L]),
                  uplo = "U", diag = "N", x = lu.A@x))
    ae1(A, with(e.lu.A, P1. \%*\% L \%*\% U))
    ae2(A[asPerm(lu.A@perm), ], with(e.lu.A, L \%*\% U))
})

## Factorization handled as factorized matrix
b <- rnorm(n)
stopifnot(identical(det(A), det(lu.A)),
          identical(solve(A, b), solve(lu.A, b)))
}
