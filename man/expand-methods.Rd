\name{expand-methods}
\title{Expand Matrix Factorizations}
%
\docType{methods}
\keyword{algebra}
\keyword{array}
\keyword{methods}
%
\alias{expand}
\alias{expand-methods}
\alias{expand1}
\alias{expand1-methods}
\alias{expand2}
\alias{expand2-methods}
%
\alias{expand,CHMfactor-method}
\alias{expand,denseLU-method}
\alias{expand,sparseLU-method}
%
\alias{expand1,BunchKaufman-method}
\alias{expand1,CHMsimpl-method}
\alias{expand1,CHMsuper-method}
\alias{expand1,Cholesky-method}
\alias{expand1,Schur-method}
\alias{expand1,denseLU-method}
\alias{expand1,pBunchKaufman-method}
\alias{expand1,pCholesky-method}
\alias{expand1,sparseLU-method}
\alias{expand1,sparseQR-method}
%
\alias{expand2,BunchKaufman-method}
\alias{expand2,CHMsimpl-method}
\alias{expand2,CHMsuper-method}
\alias{expand2,Cholesky-method}
\alias{expand2,Schur-method}
\alias{expand2,denseLU-method}
\alias{expand2,pBunchKaufman-method}
\alias{expand2,pCholesky-method}
\alias{expand2,sparseLU-method}
\alias{expand2,sparseQR-method}
%
\description{
  \code{expand1} and \code{expand2} construct matrix factors from
  objects specifying matrix factorizations.  Such objects typically
  do not store the factors explicitly, employing instead a compact
  representation to save memory.
}
\usage{
expand1(x, which, \dots)
expand2(x, \dots)

expand (x, \dots)
}
\arguments{
  \item{x}{a matrix factorization, typically inheriting from
    virtual class \code{\linkS4class{MatrixFactorization}}.}
  \item{which}{a character string indicating a matrix factor.}
  \item{\dots}{further arguments passed to or from methods.}
}
\value{
  \code{expand1} returns an object inheriting from virtual class
  \code{\linkS4class{Matrix}}, representing the factor indicated
  by \code{which}, always without row and column names.

  \code{expand2} returns a list of factors, typically with names
  using conventional notation, as in \code{list(L=, U=)}.
  The first and last factors get the row and column names of the
  factorized matrix, which are preserved in the \code{Dimnames}
  slot of \code{x}.
}
\details{
  Methods for \code{expand} are retained only for backwards
  compatibility with \pkg{Matrix} \code{< 1.6-0}.  New code
  should use \code{expand1} and \code{expand2}, whose methods
  provide more control and behave more consistently.  Notably,
  \code{expand2} obeys the rule that the product of the matrix
  factors in the returned list should reproduce
  (within some tolerance) the factorized matrix,
  \emph{including} its \code{dimnames}.
  
  Hence if \code{x} is a matrix and \code{y} is its factorization,
  then
  \preformatted{    all.equal(as(x, "matrix"), as(Reduce(`\%*\%`, expand2(y)), "matrix"))}
  should in most cases return \code{TRUE}.
}
\section{Methods}{
  The following table lists methods for \code{expand1} together with
  allowed values of argument \code{which}.

  \tabular{rl}{
    \code{class(x)} \tab \code{which}\cr
    \code{\linkS4class{Schur}} \tab \code{c("Q", "T", "Q.")}\cr
    \code{\linkS4class{denseLU}} \tab \code{c("P1", "P1.", "L", "U")}\cr
    \code{\linkS4class{sparseLU}} \tab \code{c("P1", "P1.", "P2", "P2.", "L", "U")}\cr
    \code{\linkS4class{sparseQR}} \tab \code{c("P1", "P1.", "P2", "P2.", "Q", "Q1", "R", "R1")}\cr
    \code{\linkS4class{BunchKaufman}}, \code{\linkS4class{pBunchKaufman}} \tab \code{c("U", "DU", "U.", "L", "DL", "L.")}\cr
    \code{\linkS4class{Cholesky}}, \code{\linkS4class{pCholesky}} \tab \code{c("P1", "P1.", "L1", "D", "L1.", "L", "L.")}\cr
    \code{\linkS4class{CHMsimpl}}, \code{\linkS4class{CHMsimpl}} \tab \code{c("P1", "P1.", "L1", "D", "L1.", "L", "L.")}
  }

  Methods for \code{expand2} and \code{expand} are described
  below.  Factor names and classes apply also to \code{expand1}.
  
  \describe{
    \item{\code{expand2}}{\code{signature(x = "CHMsimpl")}:
      expands the factorization
      \eqn{A = P_{1}' L_{1} D L_{1}' P_{1} = P_{1}' L L' P_{1}}{A = P1' * L1 * D * L1' * P1 = P1' * L * L' * P1}
      as \code{list(P1., L1, D, L1., P1)} (the default)
      or as \code{list(P1., L, L., P1)},
      depending on optional logical argument \code{LDL}.
      \code{P1} and \code{P1.} are \code{\linkS4class{pMatrix}},
      \code{L1}, \code{L1.}, \code{L}, and \code{L.} are
      \code{\linkS4class{dtCMatrix}},
      and \code{D} is a \code{\linkS4class{ddiMatrix}}.}
    \item{\code{expand2}}{\code{signature(x = "CHMsuper")}:
      as \code{CHMsimpl}, but the triangular factors are
      stored as \code{\linkS4class{dgCMatrix}}.}
    \item{\code{expand2}}{\code{signature(x = "p?Cholesky")}:
      expands the factorization
      \eqn{A = L_{1} D L_{1}' = L L'}{A = L1 * D * L1' = L * L'}
      as \code{list(L1, D, L1.)} (the default) or as \code{list(L, L.)},
      depending on optional logical argument \code{LDL}.
      \code{L1}, \code{L1.}, \code{L}, and \code{L.} are
      \code{\linkS4class{dtrMatrix}} or \code{\linkS4class{dtpMatrix}},
      and \code{D} is a \code{\linkS4class{ddiMatrix}}.}
    \item{\code{expand2}}{\code{signature(x = "p?BunchKaufman")}:
      expands the factorization
      \eqn{A = U D_{U} U' = L D_{L} L'}{A = U * DU * U' = L * DL * L'}
      where
      \eqn{U = \prod_{k = 1}^{b_{U}} P_{k} U_{k}}{U = prod(Pk * Uk : k = 1,...,bU)}
      and
      \eqn{L = \prod_{k = 1}^{b_{L}} P_{k} L_{k}}{L = prod(Pk * Lk : k = 1,...,bL)}
      as \code{list(U, DU, U.)} or \code{list(L, DL, L.)},
      depending on \code{x@uplo}.  If optional argument \code{complete}
      is \code{TRUE}, then an unnamed list giving the full expansion
      with \eqn{2 b_{U} + 1}{2*bU+1} or \eqn{2 b_{L} + 1}{2*bL+1} matrix
      factors is returned instead.
      \eqn{P_{k}}{Pk} are represented as \code{\linkS4class{pMatrix}},
      \eqn{U_{k}}{Uk} and \eqn{L_{k}}{Lk} are represented as
      \code{\linkS4class{dtCMatrix}}, and
      \eqn{D_{U}}{DU} and \eqn{D_{L}}{DL} are represented as
      \code{\linkS4class{dsCMatrix}}.}
    \item{\code{expand2}}{\code{signature(x = "Schur")}:
      expands the factorization
      \eqn{A = Q T Q'}{A = Q * T * Q'}
      as \code{list(Q, T, Q.)}.
      \code{Q} and \code{Q.} are \code{x@Q} and \code{t(x@Q)}
      modulo \code{Dimnames}, and \code{T} is \code{x@T}.}
    \item{\code{expand2}}{\code{signature(x = "sparseLU")}:
      expands the factorization
      \eqn{A = P_{1}' L U P_{2}'}{A = P1' * L * U * P2'}
      as \code{list(P1., L, U, P2.)}.
      \code{P1.} and \code{P2.} are \code{\linkS4class{pMatrix}},
      and \code{L} and \code{U} are \code{\linkS4class{dtCMatrix}}.}
    \item{\code{expand2}}{\code{signature(x = "denseLU")}:
      expands the factorization
      \eqn{A = P_{1}' L U}{A = P1' * L * U}
      as \code{list(P1., L, U)}.
      \code{P1.} is a \code{\linkS4class{pMatrix}},
      and \code{L} and \code{U} are \code{\linkS4class{dtrMatrix}}
      if square and \code{\linkS4class{dgeMatrix}} otherwise.}
    \item{\code{expand2}}{\code{signature(x = "sparseQR")}:
      expands the factorization
      \eqn{A = P_{1}' Q R P_{2}' = P_{1}' Q_{1} R_{1} P_{2}'}{A = P1' * Q * R * P2' = P1' * Q1 * R1 * P2'}
      as \code{list(P1., Q, R, P2.)} or \code{list(P1., Q1, R1, P2.)},
      depending on optional logical argument \code{complete}.
      \code{P1.} and \code{P2.} are \code{\linkS4class{pMatrix}},
      \code{Q} and \code{Q1} are \code{\linkS4class{dgeMatrix}},
      \code{R} is a \code{\linkS4class{dgCMatrix}},
      and \code{R1} is a \code{\linkS4class{dtCMatrix}}.}
    \item{\code{expand}}{\code{signature(x = "CHMfactor")}:
      as \code{expand2}, but returning \code{list(P, L)}.
      \code{expand(x)[["P"]]} and \code{expand2(x)[["P1"]]}
      represent the same permutation matrix \eqn{P_{1}}{P1}
      but have opposite \code{margin} slots and inverted
      \code{perm} slots.  The components of \code{expand(x)}
      do not preserve \code{x@Dimnames}.}
    \item{\code{expand}}{\code{signature(x = "sparseLU")}:
      as \code{expand2}, but returning \code{list(P, L, U, Q)}.
      \code{expand(x)[["Q"]]} and \code{expand2(x)[["P2."]]}
      represent the same permutation matrix \eqn{P_{2}'}{P2'}
      but have opposite \code{margin} slots and inverted
      \code{perm} slots.  \code{expand(x)[["P"]]} represents
      the permutation matrix \eqn{P_{1}}{P1} rather than its
      transpose \eqn{P_{1}'}{P1'}; it is \code{expand2(x)[["P1."]]}
      with an inverted \code{perm} slot.  \code{expand(x)[["L"]]}
      and \code{expand2(x)[["L"]]} represent the same unit lower
      triangular matrix \eqn{L}, but with \code{diag} slot equal
      to \code{"N"} and \code{"U"}, respectively.
      \code{expand(x)[["L"]]} and \code{expand(x)[["U"]]}
      store the \emph{permuted} first and second components of
      \code{x@Dimnames} in their \code{Dimnames} slots.}
    \item{\code{expand}}{\code{signature(x = "denseLU")}:
      as \code{expand2}, but returning \code{list(L, U, P)}.
      \code{expand(x)[["P"]]} and \code{expand2(x)[["P1."]]}
      are identical modulo \code{Dimnames}.  The components
      of \code{expand(x)} do not preserve \code{x@Dimnames}.}
  }
}
\seealso{
  The virtual class \code{\linkS4class{MatrixFactorization}}
  of matrix factorizations.
  
  Generic functions \code{\link{Cholesky}}, \code{\link{BunchKaufman}},
  \code{\link{Schur}}, \code{\link{lu}}, and \code{\link{qr}} for
  \emph{computing} factorizations.
}
\examples{
\dontshow{ % for R_DEFAULT_PACKAGES=NULL
library(stats, pos = "package:base", verbose = FALSE)
}
showMethods("expand1", inherited = FALSE)
showMethods("expand2", inherited = FALSE)
set.seed(0)

(A <- Matrix(rnorm(9L, 0, 10), 3L, 3L))
(lu.A <- lu(A))
(e.lu.A <- expand2(lu.A))
stopifnot(exprs = {
    is.list(e.lu.A)
    identical(names(e.lu.A), c("P1.", "L", "U"))
    all(sapply(e.lu.A, is, "Matrix"))
    all.equal(as(A, "matrix"), as(Reduce(`\%*\%`, e.lu.A), "matrix"))
})

## 'expand1' and 'expand2' give equivalent results modulo
## dimnames and representation of permutation matrices;
## see also function 'alt' in example("Cholesky-methods")
(a1 <- sapply(names(e.lu.A), expand1, x = lu.A, simplify = FALSE))
all.equal(a1, e.lu.A)

## see help("denseLU-class") and others for more examples
}
