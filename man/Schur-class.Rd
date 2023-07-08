\name{Schur-class}
\title{Schur Factorizations}
%
\docType{class}
\keyword{algebra}
\keyword{array}
\keyword{classes}
%
\alias{Schur-class}
%
\alias{determinant,Schur,logical-method}
%
\description{
  \code{Schur} is the class of Schur factorizations of
  \eqn{n \times n}{n-by-n} real matrices \eqn{A},
  having the general form
  \deqn{A = Q T Q'}{A = Q * T * Q'}
  where
  \eqn{Q} is an orthogonal matrix and
  \eqn{T} is a block upper triangular matrix with
  \eqn{1 \times 1}{1-by-1} or \eqn{2 \times 2}{2-by-2} diagonal blocks
  specifying the real and complex conjugate eigenvalues of \eqn{A}.
  The column vectors of \eqn{Q} are the Schur vectors of \eqn{A},
  and \eqn{T} is the Schur form of \eqn{A}.

  The Schur factorization generalizes the spectral decomposition
  of normal matrices \eqn{A}, whose Schur form is block diagonal,
  to arbitrary square matrices.
}
\section{Slots}{
  \describe{
    \item{\code{Dim}, \code{Dimnames}}{inherited from virtual class
      \code{\linkS4class{MatrixFactorization}}.}
    \item{\code{Q}}{an orthogonal matrix,
      inheriting from virtual class \code{\linkS4class{Matrix}}.}
    \item{\code{T}}{a block upper triangular matrix,
      inheriting from virtual class \code{\linkS4class{Matrix}}.
      The diagonal blocks have dimensions 1-by-1 or 2-by-2.}
    \item{\code{EValues}}{a numeric or complex vector containing
      the eigenvalues of the diagonal blocks of \code{T}, which are
      the eigenvalues of \code{T} and consequently of the factorized
      matrix.}
  }
}
\section{Extends}{
  Class \code{\linkS4class{SchurFactorization}}, directly.
  Class \code{\linkS4class{MatrixFactorization}}, by class
  \code{\linkS4class{SchurFactorization}}, distance 2.
}
\section{Instantiation}{
  Objects can be generated directly by calls of the form
  \code{new("Schur", ...)}, but they are more typically obtained
  as the value of \code{\link{Schur}(x)} for \code{x} inheriting from
  \code{\linkS4class{Matrix}} (often \code{\linkS4class{dgeMatrix}}).
}
\section{Methods}{
  \describe{
    \item{\code{determinant}}{\code{signature(from = "Schur", logarithm = "logical")}:
      computes the determinant of the factorized matrix \eqn{A}
      or its logarithm.}
    \item{\code{expand1}}{\code{signature(x = "Schur")}:
      see \code{\link{expand1-methods}}.}
    \item{\code{expand2}}{\code{signature(x = "Schur")}:
      see \code{\link{expand2-methods}}.}
    \item{\code{solve}}{\code{signature(a = "Schur", b = .)}:
      see \code{\link{solve-methods}}.}
  }
}
\details{
  The matrix \eqn{A} and its Schur form \eqn{T} are \emph{similar}
  and thus have the same spectrum.  The eigenvalues are computed
  trivially as the eigenvalues of the diagonal blocks of \eqn{T}.
}
\seealso{
  Class \code{\linkS4class{dgeMatrix}}.
  
  Generic functions \code{\link{Schur}},
  \code{\link{expand1}} and \code{\link{expand2}}.
}
\references{
  The LAPACK source code, including documentation; see
  \url{https://netlib.org/lapack/double/dgees.f}.

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
showClass("Schur")
set.seed(0)

n <- 4L
(A <- Matrix(rnorm(n * n), n, n))

## With dimnames, to see that they are propagated :
dimnames(A) <- list(paste0("r", seq_len(n)),
                    paste0("c", seq_len(n)))

(sch.A <- Schur(A))
str(e.sch.A <- expand2(sch.A), max.level = 2L)

## A ~ Q T Q' in floating point
stopifnot(exprs = {
    identical(names(e.sch.A), c("Q", "T", "Q."))
    all.equal(A, with(e.sch.A, Q \%*\% T \%*\% Q.))
})

## Factorization handled as factorized matrix
b <- rnorm(n)
stopifnot(all.equal(det(A), det(sch.A)),
          all.equal(solve(A, b), solve(sch.A, b)))

## One of the non-general cases:
Schur(Diagonal(6L))
}
