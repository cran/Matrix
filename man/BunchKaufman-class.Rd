\name{BunchKaufman-class}
\title{Dense Bunch-Kaufman Factorizations}
%
\docType{class}
\keyword{algebra}
\keyword{array}
\keyword{classes}
%
\alias{BunchKaufman-class}
\alias{pBunchKaufman-class}
%
\alias{coerce,BunchKaufman,dtrMatrix-method}
\alias{coerce,pBunchKaufman,dtpMatrix-method}
\alias{determinant,BunchKaufman,logical-method}
\alias{determinant,pBunchKaufman,logical-method}
%
\description{
  Classes \code{BunchKaufman} and \code{pBunchKaufman} represent
  Bunch-Kaufman factorizations of \eqn{n \times n}{n-by-n} real,
  symmetric matrices \eqn{A}, having the general form
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
  
  These classes store the nonzero entries of the
  \eqn{2 b_{U} + 1}{2*bU+1} or \eqn{2 b_{L} + 1}{2*bL+1} factors,
  which are individually sparse,
  in a dense format as a vector of length
  \eqn{nn}{n*n} (\code{BunchKaufman}) or
  \eqn{n(n+1)/2}{n*(n+1)/2} (\code{pBunchKaufman}),
  the latter giving the \dQuote{packed} representation.
}
\section{Slots}{
  \describe{
    \item{\code{Dim}, \code{Dimnames}}{inherited from virtual class
      \code{\linkS4class{MatrixFactorization}}.}
    \item{\code{uplo}}{a string, either \code{"U"} or \code{"L"},
      indicating which triangle (upper or lower) of the factorized
      symmetric matrix was used to compute the factorization and
      in turn how the \code{x} slot is partitioned.}
    \item{\code{x}}{a numeric vector of length \code{n*n}
      (\code{BunchKaufman}) or \code{n*(n+1)/2} (\code{pBunchKaufman}),
      where \code{n=Dim[1]}.
      The details of the representation are specified by the manual
      for LAPACK routines \code{dsytrf} and \code{dsptrf}.}
    \item{\code{perm}}{an integer vector of length \code{n=Dim[1]}
      specifying row and column interchanges as described in the manual
      for LAPACK routines \code{dsytrf} and \code{dsptrf}.}
  }
}
\section{Extends}{
  Class \code{\linkS4class{BunchKaufmanFactorization}}, directly.
  Class \code{\linkS4class{MatrixFactorization}}, by class
  \code{\linkS4class{BunchKaufmanFactorization}}, distance 2.
}
\section{Instantiation}{
  Objects can be generated directly by calls of the form
  \code{new("BunchKaufman", ...)} or \code{new("pBunchKaufman", ...)},
  but they are more typically obtained as the value of
  \code{\link{BunchKaufman}(x)} for \code{x} inheriting from
  \code{\linkS4class{dsyMatrix}} or \code{\linkS4class{dspMatrix}}.
}
\section{Methods}{
  \describe{
    \item{\code{coerce}}{\code{signature(from = "BunchKaufman", to = "dtrMatrix")}:
      returns a \code{\linkS4class{dtrMatrix}}, useful for inspecting
      the internal representation of the factorization; see \sQuote{Note}.}
    \item{\code{coerce}}{\code{signature(from = "pBunchKaufman", to = "dtpMatrix")}:
      returns a \code{\linkS4class{dtpMatrix}}, useful for inspecting
      the internal representation of the factorization; see \sQuote{Note}.}
    \item{\code{determinant}}{\code{signature(from = "p?BunchKaufman", logarithm = "logical")}:
      computes the determinant of the factorized matrix \eqn{A}
      or its logarithm.}
    \item{\code{expand1}}{\code{signature(x = "p?BunchKaufman")}:
      see \code{\link{expand1-methods}}.}
    \item{\code{expand2}}{\code{signature(x = "p?BunchKaufman")}:
      see \code{\link{expand2-methods}}.}
    \item{\code{solve}}{\code{signature(a = "p?BunchKaufman", b = .)}:
      see \code{\link{solve-methods}}.}
  }
}
\note{
  In \pkg{Matrix} \code{< 1.6-0}, class \code{BunchKaufman} extended
  \code{\linkS4class{dtrMatrix}} and class \code{pBunchKaufman} extended
  \code{\linkS4class{dtpMatrix}}, reflecting the fact that the internal
  representation of the factorization is fundamentally triangular:
  there are \eqn{n(n+1)/2}{n*(n+1)/2} \dQuote{parameters}, and these
  can be arranged systematically to form an \eqn{n \times n}{n-by-n}
  triangular matrix.
  \pkg{Matrix} \code{1.6-0} removed these extensions so that methods
  would no longer be inherited from \code{dtrMatrix} and \code{dtpMatrix}.  
  The availability of such methods gave the wrong impression that
  \code{BunchKaufman} and \code{pBunchKaufman} represent a (singular)
  matrix, when in fact they represent an ordered set of matrix factors.
  
  The coercions \code{as(., "dtrMatrix")} and \code{as(., "dtpMatrix")}
  are provided for users who understand the caveats.
}
\seealso{
  Class \code{\linkS4class{dsyMatrix}} and its packed counterpart.

  Generic functions \code{\link{BunchKaufman}},
  \code{\link{expand1}}, and \code{\link{expand2}}.
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
library(stats, pos = "package:base", verbose = FALSE)
library(utils, pos = "package:base", verbose = FALSE)
}
showClass("BunchKaufman")
set.seed(1)

n <- 6L
(A <- forceSymmetric(Matrix(rnorm(n * n), n, n)))

## With dimnames, to see that they are propagated :
dimnames(A) <- rep.int(list(paste0("x", seq_len(n))), 2L)

(bk.A <- BunchKaufman(A))
str(e.bk.A <- expand2(bk.A, complete = FALSE), max.level = 2L)
str(E.bk.A <- expand2(bk.A, complete =  TRUE), max.level = 2L)

## Underlying LAPACK representation
(m.bk.A <- as(bk.A, "dtrMatrix"))
stopifnot(identical(as(m.bk.A, "matrix"), `dim<-`(bk.A@x, bk.A@Dim)))

## Number of factors is 2*b+1, b <= n, which can be nontrivial ...
(b <- (length(E.bk.A) - 1L) \%/\% 2L)

ae1 <- function(a, b, ...) all.equal(as(a, "matrix"), as(b, "matrix"), ...)
ae2 <- function(a, b, ...) ae1(unname(a), unname(b), ...)

## A ~ U DU U', U := prod(Pk Uk) in floating point
stopifnot(exprs = {
    identical(names(e.bk.A), c("U", "DU", "U."))
    identical(e.bk.A[["U" ]], Reduce(`\%*\%`, E.bk.A[seq_len(b)]))
    identical(e.bk.A[["U."]], t(e.bk.A[["U"]]))
    ae1(A, with(e.bk.A, U \%*\% DU \%*\% U.))
})

## Factorization handled as factorized matrix
b <- rnorm(n)
stopifnot(identical(det(A), det(bk.A)),
          identical(solve(A, b), solve(bk.A, b)))
}
