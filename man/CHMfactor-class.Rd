\name{CHMfactor-class}
\title{Sparse Cholesky Factorizations}
%
\docType{class}
\keyword{algebra}
\keyword{array}
\keyword{classes}
\keyword{programming}
\keyword{utilities}
%
\alias{CHMfactor-class}
\alias{CHMsimpl-class}
\alias{CHMsuper-class}
\alias{dCHMsimpl-class}
\alias{dCHMsuper-class}
\alias{nCHMsimpl-class}
\alias{nCHMsuper-class}
%
\alias{coerce,CHMsimpl,dtCMatrix-method}
\alias{coerce,CHMsuper,dgCMatrix-method}
\alias{determinant,CHMfactor,logical-method}
\alias{diag,CHMfactor-method}
\alias{update,CHMfactor-method}
%
\alias{isLDL}
%
\description{
  \code{CHMfactor} is the virtual class of sparse Cholesky
  factorizations of \eqn{n \times n}{n-by-n} real, symmetric
  matrices \eqn{A}, having the general form
  \deqn{P_1 A P_1' = L_1 D L_1' \overset{D_{jj} \ge 0}{=} L L'}{P1 * A * P1' = L1 * D * L1' [ = L * L' ]}
  or (equivalently)
  \deqn{A = P_1' L_1 D L_1' P_1 \overset{D_{jj} \ge 0}{=} P_1' L L' P_1}{A = P1' L1 * D * L1' * P1 [ = P1' * L * L' * P1 ]}
  where
  \eqn{P_1}{P1} is a permutation matrix,
  \eqn{L_1}{L1} is a unit lower triangular matrix,
  \eqn{D} is a diagonal matrix, and
  \eqn{L = L_1 \sqrt{D}}{L = L1 * sqrt(D)}.
  The second equalities hold only for positive semidefinite \eqn{A},
  for which the diagonal entries of \eqn{D} are non-negative
  and \eqn{\sqrt{D}}{sqrt(D)} is well-defined.
  
  The implementation of class \code{CHMfactor} is based on
  CHOLMOD's C-level \code{cholmod_factor_struct}.  Virtual
  subclasses \code{CHMsimpl} and \code{CHMsuper} separate
  the simplicial and supernodal variants.  These have nonvirtual
  subclasses \code{[dn]CHMsimpl} and \code{[dn]CHMsuper},
  where prefix \samp{d} and prefix \samp{n} are reserved
  for numeric and symbolic factorizations, respectively.
}
\usage{
isLDL(x)
}
\arguments{
  \item{x}{an object inheriting from virtual class \code{CHMfactor},
    almost always the result of a call to generic function
    \code{\link{Cholesky}}.}
}
\value{
  \code{isLDL(x)} returns \code{TRUE} or \code{FALSE}:
  \code{TRUE} if \code{x} stores the lower triangular entries
  of \eqn{L_1-I+D}{L1-I+D},
  \code{FALSE} if \code{x} stores the lower triangular entries
  of \eqn{L}.
}
\section{Slots}{
  Of \code{CHMfactor}:
  \describe{
    \item{\code{Dim}, \code{Dimnames}}{inherited from virtual class
      \code{\linkS4class{MatrixFactorization}}.}
    \item{\code{colcount}}{an integer vector of length \code{Dim[1]}
      giving an \emph{estimate} of the number of nonzero entries in
      each column of the lower triangular Cholesky factor.
      If symbolic analysis was performed prior to factorization,
      then the estimate is exact.}
    \item{\code{perm}}{a 0-based integer vector of length \code{Dim[1]}
      specifying the permutation applied to the rows and columns
      of the factorized matrix.  \code{perm} of length 0 is valid and
      equivalent to the identity permutation, implying no pivoting.}
    \item{\code{type}}{an integer vector of length 6 specifying
      details of the factorization.  The elements correspond to
      members \code{ordering}, \code{is_ll}, \code{is_super},
      \code{is_monotonic}, \code{maxcsize}, and \code{maxesize}
      of the original \code{cholmod_factor_struct}.
      Simplicial and supernodal factorizations are distinguished
      by \code{is_super}.  Simplicial factorizations do not use
      \code{maxcsize} or \code{maxesize}.  Supernodal factorizations
      do not use \code{is_ll} or \code{is_monotonic}.}
  }
  Of \code{CHMsimpl} (all unused by \code{nCHMsimpl}): % FIXME
  \describe{
    \item{\code{nz}}{an integer vector of length \code{Dim[1]}
      giving the number of nonzero entries in each column of the
      lower triangular Cholesky factor.  There is at least one
      nonzero entry in each column, because the diagonal elements
      of the factor are stored explicitly.}
    \item{\code{p}}{an integer vector of length \code{Dim[1]+1}.
      Row indices of nonzero entries in column \code{j} of the
      lower triangular Cholesky factor are obtained as
      \code{i[p[j]+seq_len(nz[j])]+1}.}
    \item{\code{i}}{an integer vector of length greater than or equal
      to \code{sum(nz)} containing the row indices of nonzero entries
      in the lower triangular Cholesky factor.  These are grouped by
      column and sorted within columns, but the columns themselves
      need not be ordered monotonically.  Columns may be overallocated,
      i.e., the number of elements of \code{i} reserved for column
      \code{j} may exceed \code{nz[j]}.}
    \item{\code{prv}, \code{nxt}}{integer vectors of length
      \code{Dim[1]+2} indicating the order in which the columns of
      the lower triangular Cholesky factor are stored in \code{i}
      and \code{x}.
      Starting from \code{j <- Dim[1]+2},
      the recursion \code{j <- nxt[j+1]+1} traverses the columns
      in  forward order and terminates when \code{nxt[j+1] = -1}.
      Starting from \code{j <- Dim[1]+1},
      the recursion \code{j <- prv[j+1]+1} traverses the columns
      in backward order and terminates when \code{prv[j+1] = -1}.}
  }
  Of \code{dCHMsimpl}:
  \describe{
    \item{\code{x}}{a numeric vector parallel to \code{i} containing
      the corresponding nonzero entries of the lower triangular
      Cholesky factor \eqn{L} \emph{or} (if and only if \code{type[2]}
      is 0) of the lower triangular matrix \eqn{L_1-I+D}{L1-I+D}.}
    }
  Of \code{CHMsuper}:
  \describe{
    \item{\code{super}, \code{pi}, \code{px}}{integer vectors of
      length \code{nsuper+1}, where \code{nsuper} is the number of
      supernodes.  \code{super[j]+1} is the index of the leftmost
      column of supernode \code{j}.  The row indices of supernode
      \code{j} are obtained as \code{s[pi[j]+seq_len(pi[j+1]-pi[j])]+1}.
      The numeric entries of supernode \code{j} are obtained as
      \code{x[px[j]+seq_len(px[j+1]-px[j])]+1} (if slot \code{x}
      is available).}
    \item{\code{s}}{an integer vector of length greater than or equal
      to \code{Dim[1]} containing the row indices of the supernodes.
      \code{s} may contain duplicates, but not within a supernode,
      where the row indices must be increasing.}
  }
  Of \code{dCHMsuper}:
  \describe{
    \item{\code{x}}{a numeric vector of length less than or equal to
      \code{prod(Dim)} containing the numeric entries of the supernodes.}
  }
}
\section{Extends}{
  Class \code{\linkS4class{MatrixFactorization}}, directly.
}
\section{Instantiation}{
  Objects can be generated directly by calls of the form
  \code{new("dCHMsimpl", ...)}, etc., but \code{dCHMsimpl} and
  \code{dCHMsuper} are more typically obtained as the value of
  \code{\link{Cholesky}(x, ...)} for \code{x} inheriting from
  \code{\linkS4class{sparseMatrix}}
  (often \code{\linkS4class{dsCMatrix}}).
  
  There is currently no API outside of calls to \code{\link{new}}
  for generating \code{nCHMsimpl} and \code{nCHMsuper}.  These
  classes are vestigial and may be formally deprecated in a future
  version of \pkg{Matrix}.
}
\section{Methods}{
  \describe{
    \item{\code{coerce}}{\code{signature(from = "CHMsimpl", to = "dtCMatrix")}:
      returns a \code{\linkS4class{dtCMatrix}} representing
      the lower triangular Cholesky factor \eqn{L} \emph{or}
      the lower triangular matrix \eqn{L_1-I+D}{L1-I+D},
      the latter if and only if \code{from@type[2]} is 0.}
    \item{\code{coerce}}{\code{signature(from = "CHMsuper", to = "dgCMatrix")}:
      returns a \code{\linkS4class{dgCMatrix}} representing
      the lower triangular Cholesky factor \eqn{L}.  Note that,
      for supernodes spanning two or more columns, the supernodal
      algorithm by design stores non-structural zeros above
      the main diagonal, hence \code{\linkS4class{dgCMatrix}} is
      indeed more appropriate than \code{\linkS4class{dtCMatrix}}
      as a coercion target.}
    \item{\code{determinant}}{\code{signature(from = "CHMfactor", logarithm = "logical")}:
      behaves according to an optional argument \code{sqrt}.
      If \code{sqrt = FALSE}, then this method computes the determinant
      of the factorized matrix \eqn{A} or its logarithm.
      If \code{sqrt = TRUE}, then this method computes the determinant
      of the factor \eqn{L = L_1 sqrt(D)}{L = L1 * sqrt(D)} or
      its logarithm, giving \code{NaN} for the modulus when \eqn{D}
      has negative diagonal elements.  For backwards compatibility,
      the default value of \code{sqrt} is \code{TRUE}, but that can
      be expected change in a future version of \pkg{Matrix}, hence
      defensive code will always set \code{sqrt} (to \code{TRUE},
      if the code must remain backwards compatible with \pkg{Matrix}
      \code{< 1.6-0}).  Calls to this method not setting \code{sqrt}
      may warn about the pending change.  The warnings can be disabled
      with \code{\link{options}(Matrix.warnSqrtDefault = 0)}.}
    \item{\code{diag}}{\code{signature(x = "CHMfactor")}:
      returns a numeric vector of length \eqn{n} containing the diagonal
      elements of \eqn{D}, which (\emph{if} they are all non-negative)
      are the squared diagonal elements of \eqn{L}.}
    \item{\code{expand}}{\code{signature(x = "CHMfactor")}:
      see \code{\link{expand-methods}}.}
    \item{\code{expand1}}{\code{signature(x = "CHMsimpl")}:
      see \code{\link{expand1-methods}}.}
    \item{\code{expand1}}{\code{signature(x = "CHMsuper")}:
      see \code{\link{expand1-methods}}.}
    \item{\code{expand2}}{\code{signature(x = "CHMsimpl")}:
      see \code{\link{expand2-methods}}.}
    \item{\code{expand2}}{\code{signature(x = "CHMsuper")}:
      see \code{\link{expand2-methods}}.}
    \item{\code{image}}{\code{signature(x = "CHMfactor")}:
      see \code{\link{image-methods}}.}
    \item{\code{nnzero}}{\code{signature(x = "CHMfactor")}:
      see \code{\link{nnzero-methods}}.}
    \item{\code{solve}}{\code{signature(a = "CHMfactor", b = .)}:
      see \code{\link{solve-methods}}.}
    \item{\code{update}}{\code{signature(object = "CHMfactor")}:
      returns a copy of \code{object} with the same nonzero pattern
      but with numeric entries updated according to additional
      arguments \code{parent} and \code{mult}, where \code{parent}
      is (coercible to) a \code{\linkS4class{dsCMatrix}} or a
      \code{\linkS4class{dgCMatrix}} and \code{mult} is a numeric
      vector of positive length.
      \cr
      The numeric entries are updated with those of the Cholesky
      factor of \code{F(parent) + mult[1] * I}, i.e.,
      \code{F(parent)} plus \code{mult[1]} times the identity matrix,
      where \code{F = \link{identity}} for symmetric \code{parent}
      and \code{F = \link{tcrossprod}} for other \code{parent}.
      The nonzero pattern of \code{F(parent)} must match
      that of \code{S} if \code{object = Cholesky(S, ...)}.}
    \item{\code{updown}}{\code{signature(update = ., C = ., object = "CHMfactor")}:
      see \code{\link{updown-methods}}.}
  }
}
\seealso{
  Class \code{\linkS4class{dsCMatrix}}.

  Generic functions \code{\link{Cholesky}}, \code{\link{updown}},
  \code{\link{expand1}} and \code{\link{expand2}}.
}
\references{
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
library(utils, pos = "package:base", verbose = FALSE)
}
showClass("dCHMsimpl")
showClass("dCHMsuper")
set.seed(2)

m <- 1000L
n <- 200L
M <- rsparsematrix(m, n, 0.01)
A <- crossprod(M)

## With dimnames, to see that they are propagated :
dimnames(A) <- dn <- rep.int(list(paste0("x", seq_len(n))), 2L)

(ch.A <- Cholesky(A)) # pivoted, by default
str(e.ch.A <- expand2(ch.A, LDL =  TRUE), max.level = 2L)
str(E.ch.A <- expand2(ch.A, LDL = FALSE), max.level = 2L)

ae1 <- function(a, b, ...) all.equal(as(a, "matrix"), as(b, "matrix"), ...)
ae2 <- function(a, b, ...) ae1(unname(a), unname(b), ...)

## A ~ P1' L1 D L1' P1 ~ P1' L L' P1 in floating point
stopifnot(exprs = {
    identical(names(e.ch.A), c("P1.", "L1", "D", "L1.", "P1"))
    identical(names(E.ch.A), c("P1.", "L" ,      "L." , "P1"))
    identical(e.ch.A[["P1"]],
              new("pMatrix", Dim = c(n, n), Dimnames = c(list(NULL), dn[2L]),
                  margin = 2L, perm = invertPerm(ch.A@perm, 0L, 1L)))
    identical(e.ch.A[["P1."]], t(e.ch.A[["P1"]]))
    identical(e.ch.A[["L1."]], t(e.ch.A[["L1"]]))
    identical(E.ch.A[["L." ]], t(E.ch.A[["L" ]]))
    identical(e.ch.A[["D"]], Diagonal(x = diag(ch.A)))
    all.equal(E.ch.A[["L"]], with(e.ch.A, L1 \%*\% sqrt(D)))
    ae1(A, with(e.ch.A, P1. \%*\% L1 \%*\% D \%*\% L1. \%*\% P1))
    ae1(A, with(E.ch.A, P1. \%*\% L  \%*\%         L.  \%*\% P1))
    ae2(A[ch.A@perm + 1L, ch.A@perm + 1L], with(e.ch.A, L1 \%*\% D \%*\% L1.))
    ae2(A[ch.A@perm + 1L, ch.A@perm + 1L], with(E.ch.A, L  \%*\%         L. ))
})

## Factorization handled as factorized matrix
## (in some cases only optionally, depending on arguments)
b <- rnorm(n)
stopifnot(identical(det(A), det(ch.A, sqrt = FALSE)),
          identical(solve(A, b), solve(ch.A, b, system = "A")))

u1 <- update(ch.A,   A , mult = sqrt(2))
u2 <- update(ch.A, t(M), mult = sqrt(2)) # updating with crossprod(M), not M
stopifnot(all.equal(u1, u2, tolerance = 1e-14))
}
