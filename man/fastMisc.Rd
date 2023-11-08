\name{fastMisc}
\title{\dQuote{Low Level} Coercions and Methods}
%
\keyword{utilities}
%
\alias{fastMisc}
% coercions:
\alias{.M2kind}
\alias{.M2gen}
\alias{.M2sym}
\alias{.M2tri}
\alias{.M2diag}
\alias{.M2v}
\alias{.M2m}
\alias{.M2unpacked}
\alias{.M2packed}
\alias{.M2C}
\alias{.M2R}
\alias{.M2T}
\alias{.M2V}
\alias{.m2V}
\alias{.sparse2dense}
\alias{.diag2dense}
\alias{.ind2dense}
\alias{.m2dense}
\alias{.dense2sparse}
\alias{.diag2sparse}
\alias{.ind2sparse}
\alias{.m2sparse}
\alias{.tCRT}
% coercions, predating API finalization, hence no longer documented:
\alias{.CR2RC}
\alias{.CR2T}
\alias{.T2CR}
\alias{.dense2g}
\alias{.dense2kind}
\alias{.dense2m}
\alias{.dense2v}
\alias{.sparse2g}
\alias{.sparse2kind}
\alias{.sparse2m}
\alias{.sparse2v}
\alias{.tCR2RC}
% other direct methods:
\alias{.diag.dsC}
\alias{.solve.dgC.lu}
\alias{.solve.dgC.qr}
\alias{.solve.dgC.chol}
\alias{.updateCHMfactor}
%
\description{
  \dQuote{Semi-API} functions used internally by \pkg{Matrix},
  often to bypass S4 dispatch and avoid the associated overhead.
  These are exported to provide this capability to expert users.
  Typical users should continue to rely on S4 generic functions
  to dispatch suitable methods, by calling,
  e.g., \code{as(., <class>)} for coercions.
}
\usage{
.M2kind(from, kind = ".", sparse = NA)

.M2gen(from, kind = ".")
.M2sym(from, \dots)
.M2tri(from, \dots)
.M2diag(from)

.M2v(from)
.M2m(from)
.M2unpacked(from)
.M2packed(from)
.M2C(from)
.M2R(from)
.M2T(from)

.M2V(from)
.m2V(from, kind = ".")

.sparse2dense(from, packed = FALSE)
.diag2dense(from, kind = ".", shape = "t", packed = FALSE, uplo = "U")
.ind2dense(from, kind = "n")
.m2dense(from, class = ".ge", uplo = "U", diag = "N", trans = FALSE)

.dense2sparse(from, repr = "C")
.diag2sparse(from, kind = ".", shape = "t", repr = "C", uplo = "U")
.ind2sparse(from, kind = "n", repr = ".")
.m2sparse(from, class = ".gC", uplo = "U", diag = "N", trans = FALSE)

.tCRT(x, lazy = TRUE)

.diag.dsC(x, Chx = Cholesky(x, LDL = TRUE), res.kind = "diag")

.solve.dgC.lu  (a, b, tol = .Machine$double.eps, check = TRUE)
.solve.dgC.qr  (a, b, order = 3L, check = TRUE)
.solve.dgC.chol(a, b, check = TRUE)

.updateCHMfactor(object, parent, mult = 0)
}
\arguments{
  \item{from, x, a, b}{a \code{\linkS4class{Matrix}}, matrix, or vector.}
  \item{kind}{a string (\code{"."}, \code{","}, \code{"n"}, \code{"l"},
    or \code{"d"}) specifying the \dQuote{kind} of the result.
    \code{"."} indicates that the kind of \code{from} should be preserved.
    \code{","} is equivalent to \code{"z"} if \code{from} is complex
    and to \code{"d"} otherwise.
    \code{"n"} indicates that the result should inherit from
    \code{\linkS4class{nMatrix}} or \code{\linkS4class{nsparseVector}}
    (and so on).}
  \item{shape}{a string (\code{"."}, \code{"g"}, \code{"s"}, or
    \code{"t"}) specifying the \dQuote{shape} of the result.  \code{"."}
    indicates that the shape of \code{from} should be preserved.
    \code{"g"} indicates that the result should inherit from
    \code{\linkS4class{generalMatrix}} (and so on).}
  \item{repr}{a string (\code{"."}, \code{"C"}, \code{"R"}, or
    \code{"T"}) specifying the sparse representation of the result.
    \code{"."} is accepted only by \code{.ind2sparse} and indicates
    the most efficient representation,
    which is \code{"C"} (\code{"R"}) for \code{margin = 2} (\code{1}).
    \code{"C"} indicates that the result should inherit from
    \code{\linkS4class{CsparseMatrix}} (and so on).}
  \item{packed}{a logical indicating if the result should
    inherit from \code{\linkS4class{packedMatrix}}
    rather than from \code{\linkS4class{unpackedMatrix}}.
    It is ignored for \code{from} inheriting from
    \code{\linkS4class{generalMatrix}}.}
  \item{sparse}{a logical indicating if the result should inherit
    from \code{\linkS4class{sparseMatrix}} rather than from
    \code{\linkS4class{denseMatrix}}.  If \code{NA}, then the result
    will be formally sparse if and only if \code{from} is.}
  \item{uplo}{a string (\code{"U"} or \code{"L"}) indicating whether
    the result (if symmetric or triangular) should store the upper or
    lower triangle of \code{from}.  The elements of \code{from} in the
    opposite triangle are ignored.}
  \item{diag}{a string (\code{"N"} or \code{"U"}) indicating whether
    the result (if triangular) should be formally nonunit or unit
    triangular.  In the unit triangular case, the diagonal elements
    of \code{from} are ignored.}
  \item{trans}{a logical indicating if the result should be a 1-row
    matrix rather than a 1-column matrix where \code{from} is a vector
    but not a matrix.}
  \item{class}{a string whose first three characters specify the class
    of the result.  It should match the pattern
    \code{"^[.nld](ge|sy|tr|sp|tp)"} for \code{.m2dense} and
    \code{"^[.nld][gst][CRT]"} for \code{.m2sparse},
    where \code{"."} in the first position is equivalent to \code{"l"}
    for logical arguments and \code{"d"} for numeric arguments.}
  \item{\dots}{optional arguments passed to \code{\link{isSymmetric}}
    or \code{\link{isTriangular}}.}
  %% .tCRT :
  \item{lazy}{a logical indicating if the transpose should be
    constructed with minimal allocation, but possibly \emph{without}
    preserving representation.}
  %% .diag.dsC :
  \item{Chx}{optionally, the \code{\link{Cholesky}(x, \dots)}
    factorization of \code{x}.  If supplied, then \code{x} is unused.}
  \item{res.kind}{a string in \code{c("trace", "sumLog", "prod", "min",
      "max", "range", "diag", "diagBack")}.}
  %% .solve.dgC.* :
  \item{tol}{see \code{\link{lu-methods}}.}
  \item{order}{see \code{\link{qr-methods}}.}
  \item{check}{a logical indicating if the first argument should be
    tested for inheritance from \code{\linkS4class{dgCMatrix}} and
    coerced if necessary.   Set to \code{FALSE} for speed only if it
    is known to already inherit from \code{\linkS4class{dgCMatrix}}.}
  %% .updateCHMfactor :
  \item{object}{a Cholesky factorization inheriting from virtual class
    \code{CHMfactor}, almost always the result of a call to generic
    function \code{\link{Cholesky}}.}
  \item{parent}{an object of class \code{\linkS4class{dsCMatrix}}
    or class \code{\linkS4class{dgCMatrix}}.}
  \item{mult}{a numeric vector of postive length.  Only the first
    element is used, and that must be finite.}
}
\details{
  Functions with names of the form \code{.<A>2<B>} implement coercions
  from virtual class A to the \dQuote{nearest} non-virtual subclass of
  virtual class B, where the virtual classes are abbreviated as follows:

  \describe{
    \item{\code{M}}{\code{\linkS4class{Matrix}}}
    \item{\code{V}}{\code{\linkS4class{sparseVector}}}
    \item{\code{m}}{matrix}
    \item{\code{v}}{vector}
    \item{\code{dense}}{\code{\linkS4class{denseMatrix}}}
    \item{\code{unpacked}}{\code{\linkS4class{unpackedMatrix}}}
    \item{\code{packed}}{\code{\linkS4class{packedMatrix}}}
    \item{\code{sparse}}{%
      \code{\linkS4class{CsparseMatrix}},
      \code{\linkS4class{RsparseMatrix}}, or
      \code{\linkS4class{TsparseMatrix}}}
    \item{\code{C}}{\code{\linkS4class{CsparseMatrix}}}
    \item{\code{R}}{\code{\linkS4class{RsparseMatrix}}}
    \item{\code{T}}{\code{\linkS4class{TsparseMatrix}}}
    \item{\code{gen}}{\code{\linkS4class{generalMatrix}}}
    \item{\code{sym}}{\code{\linkS4class{symmetricMatrix}}}
    \item{\code{tri}}{\code{\linkS4class{triangularMatrix}}}
    \item{\code{diag}}{\code{\linkS4class{diagonalMatrix}}}
    \item{\code{ind}}{\code{\linkS4class{indMatrix}}}
  }

  Abbreviations should be seen as a guide, rather than as an
  exact description of behaviour.  Notably, \code{.m2dense},
  \code{.m2sparse}, and \code{.m2V} accept vectors that are
  not matrices.

  \subsection{\code{.tCRT(x)}}{
    If \code{lazy = TRUE}, then \code{.tCRT} constructs the transpose
    of \code{x} using the most efficient representation,
    which for \samp{CRT} is \samp{RCT}.  If \code{lazy = FALSE},
    then \code{.tCRT} preserves the representation of \code{x},
    behaving as the corresponding methods for generic function \code{t}.
  }
  \subsection{\code{.diag.dsC(x)}}{
    \code{.diag.dsC} computes (or uses if \code{Chx} is supplied)
    the Cholesky factorization of \code{x} as \eqn{L D L'} in order
    to calculate one of several possible statistics from the diagonal
    entries of \eqn{D}.  See \code{res.kind} under \sQuote{Arguments}.
  }
  \subsection{\code{.solve.dgC.*(a, b)}}{
    \code{.solve.dgC.lu(a, b)} needs a square matrix \code{a}.
    \code{.solve.dgC.qr(a, b)} needs a \dQuote{long} matrix \code{a},
    with \code{nrow(a) >= ncol(a)}.
    \code{.solve.dgC.chol(a, b)} needs a \dQuote{wide} matrix \code{a},
    with \code{nrow(a) <= ncol(a)}.

    All three may be used to solve sparse linear systems directly.
    Only \code{.solve.dgC.qr} and \code{.solve.dgC.chol} be used
    to solve sparse \emph{least squares} problems.
  }
  \subsection{\code{.updateCHMfactor(object, parent, mult)}}{
    \code{.updateCHMfactor} updates \code{object} with the result
    of Cholesky factorizing
    \code{F(parent) + mult[1] * diag(nrow(parent))},
    i.e., \code{F(parent)} plus \code{mult[1]} times the identity matrix,
    where \code{F = identity} if \code{parent} is a \code{dsCMatrix}
    and \code{F = tcrossprod} if \code{parent} is a \code{dgCMatrix}.
    The nonzero pattern of \code{F(parent)} must match
    that of \code{S} if \code{object = Cholesky(S, \dots)}.
  }
}
\examples{
\dontshow{ % for R_DEFAULT_PACKAGES=NULL
library(utils, pos = "package:base", verbose = FALSE)
}
D. <- diag(x = c(1, 1, 2, 3, 5, 8))
D.0 <- Diagonal(x = c(0, 0, 0, 3, 5, 8))
S. <- toeplitz(as.double(1:6))
C. <- new("dgCMatrix", Dim = c(3L, 4L),
          p = c(0L, 1L, 1L, 1L, 3L), i = c(1L, 0L, 2L), x = c(-8, 2, 3))

stopifnot(exprs = {
    identical(.M2tri (D.), as(D., "triangularMatrix"))
    identical(.M2sym (D.), as(D.,  "symmetricMatrix"))
    identical(.M2diag(D.), as(D.,   "diagonalMatrix"))
    identical(.M2kind(C., "l"),
              as(C., "lMatrix"))
    identical(.M2kind(.sparse2dense(C.), "l"),
              as(as(C., "denseMatrix"), "lMatrix"))
    identical(.diag2sparse(D.0, ".", "t", "C"),
              .dense2sparse(.diag2dense(D.0, ".", "t", TRUE), "C"))
    identical(.M2gen(.diag2dense(D.0, ".", "s", FALSE)),
              .sparse2dense(.M2gen(.diag2sparse(D.0, ".", "s", "T"))))
    identical(S.,
              .M2m(.m2sparse(S., ".sR")))
    identical(S. * lower.tri(S.) + diag(1, 6L),
              .M2m(.m2dense (S., ".tr", "L", "U")))
    identical(.M2R(C.), .M2R(.M2T(C.)))
    identical(.tCRT(C.), .M2R(t(C.)))
})

A <- tcrossprod(C.)/6 + Diagonal(3, 1/3); A[1,2] <- 3; A
stopifnot(exprs = {
    is.numeric( x. <- c(2.2, 0, -1.2) )
    all.equal(x., .solve.dgC.lu(A, c(1,0,0), check=FALSE))
    all.equal(x., .solve.dgC.qr(A, c(1,0,0), check=FALSE))
})

## Solving sparse least squares:

X <- rbind(A, Diagonal(3)) # design matrix X (for L.S.)
Xt <- t(X)                 # *transposed*  X (for L.S.)
(y <- drop(crossprod(Xt, 1:3)) + c(-1,1)/1000) # small rand.err.
str(solveCh <- .solve.dgC.chol(Xt, y, check=FALSE)) # Xt *is* dgC..
stopifnot(exprs = {
    all.equal(solveCh$coef, 1:3, tol = 1e-3)# rel.err ~ 1e-4
    all.equal(solveCh$coef, drop(solve(tcrossprod(Xt), Xt \%*\% y)))
    all.equal(solveCh$coef, .solve.dgC.qr(X, y, check=FALSE))
})
}
