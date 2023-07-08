\name{Cholesky-methods}
\title{Methods for Cholesky Factorization}
%
\docType{methods}
\keyword{algebra}
\keyword{array}
\keyword{methods}
\concept{Choleski} % alternate English spelling
%
\alias{Cholesky}
\alias{Cholesky-methods}
%
\alias{Cholesky,ddiMatrix-method}
\alias{Cholesky,diagonalMatrix-method}
\alias{Cholesky,dsCMatrix-method}
\alias{Cholesky,dsRMatrix-method}
\alias{Cholesky,dsTMatrix-method}
\alias{Cholesky,dspMatrix-method}
\alias{Cholesky,dsyMatrix-method}
\alias{Cholesky,generalMatrix-method}
\alias{Cholesky,matrix-method}
\alias{Cholesky,symmetricMatrix-method}
\alias{Cholesky,triangularMatrix-method}
%
\alias{.SuiteSparse_version}
%
\description{
  Computes the pivoted Cholesky factorization of an
  \eqn{n \times n}{n-by-n} real, symmetric matrix \eqn{A},
  which has the general form
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
  
  Methods for \code{\linkS4class{denseMatrix}} are built on
  LAPACK routines \code{dpstrf}, \code{dpotrf}, and \code{dpptrf}.
  The latter two do not permute rows or columns,
  so that \eqn{P_1}{P1} is an identity matrix.

  Methods for \code{\linkS4class{sparseMatrix}} are built on
  CHOLMOD routines \code{cholmod_analyze} and \code{cholmod_factorize_p}.
}
\usage{
Cholesky(A, \dots)
\S4method{Cholesky}{dsyMatrix}(A, perm = TRUE, tol = -1, \dots)
\S4method{Cholesky}{dspMatrix}(A, \dots)
\S4method{Cholesky}{dsCMatrix}(A, perm = TRUE, LDL = !super, super = FALSE,
    Imult = 0, \dots)
\S4method{Cholesky}{ddiMatrix}(A, \dots)
\S4method{Cholesky}{generalMatrix}(A, uplo = "U", \dots)
\S4method{Cholesky}{triangularMatrix}(A, uplo = "U", \dots)
\S4method{Cholesky}{matrix}(A, uplo = "U", \dots)
}
\arguments{
  \item{A}{a \link[=is.finite]{finite}, symmetric matrix or
    \code{\linkS4class{Matrix}} to be factorized.  If \code{A}
    is square but not symmetric, then it will be \emph{treated}
    as symmetric; see \code{uplo}.
    Methods for dense \code{A} require positive definiteness
    when \code{perm = FALSE} and positive semidefiniteness
    when \code{perm = TRUE}.
    Methods for sparse \code{A} require positive definiteness
    when \code{LDL = TRUE} and nonzero leading principal minors
    (after pivoting) when \code{LDL = FALSE}.
    Methods for sparse, \emph{diagonal} \code{A} are an exception,
    requiring positive semidefiniteness unconditionally.}
  \item{perm}{a logical indicating if the rows and columns
    of \eqn{A} should be pivoted.  Methods for sparse \code{A}
    employ the approximate minimum degree (AMD) algorithm
    in order to reduce fill-in, i.e., without regard for
    numerical stability.
    Pivoting for sparsity may introduce nonpositive leading
    principal minors, causing the factorization to fail, in
    which case it may be necessary to set \code{perm = FALSE}.}
  \item{tol}{a \link[=is.finite]{finite} numeric tolerance,
    used only if \code{perm = TRUE}.
    The factorization algorithm stops if the pivot is less than
    or equal to \code{tol}.  Negative \code{tol} is equivalent
    to \code{nrow(A) * .Machine$double.eps * max(diag(A))}.}
  \item{LDL}{a logical indicating if the simplicial factorization
    should be computed as
    \eqn{P_1' L_1 D L_1' P_1}{P1' * L1 * D * L1' * P1},
    such that the result stores the lower triangular entries
    of \eqn{L_1-I+D}{L1-I+D}.
    The alternative is \eqn{P_1' L L' P_1}{P1' * L * L' * P1},
    such that the result stores the lower triangular entries
    of \eqn{L = L_1 \sqrt{D}}{L = L1 * sqrt(D)}.  This argument
    is ignored if \code{super = TRUE} (or if \code{super = NA}
    and the supernodal algorithm is chosen), as the supernodal
    code does not yet support the \code{LDL = TRUE} variant.}
  \item{super}{a logical indicating if the factorization should
    use the supernodal algorithm.  The alternative is the simplicial
    algorithm.  Setting \code{super = NA} leaves the choice to
    a CHOLMOD-internal heuristic.}
  \item{Imult}{a \link[=is.finite]{finite} number. The matrix
    that is factorized is \code{A + Imult * diag(nrow(A))},
    i.e., \code{A} plus \code{Imult} times the identity matrix.
    This argument is useful for symmetric, indefinite \code{A},
    as \code{Imult > max(rowSums(abs(A)) - diag(abs(A)))} ensures
    that \code{A + Imult * diag(nrow(A))} is diagonally dominant.
    (Symmetric, diagonally dominant matrices are positive definite.)}
  \item{uplo}{a string, either \code{"U"} or \code{"L"},
    indicating which triangle of \code{A} should be used
    to compute the factorization.  The default is \code{"U"},
    even for lower triangular \code{A}, to be consistent with
    \code{\link[base]{chol}} from \pkg{base}.}
  \item{\dots}{further arguments passed to or from methods.}
}
\value{
  An object representing the factorization, inheriting from
  virtual class \code{\linkS4class{CholeskyFactorization}}.
  For a traditional matrix \code{A}, the specific class is
  \code{\linkS4class{Cholesky}}.
  For \code{A} inheriting from
  \code{\linkS4class{unpackedMatrix}},
  \code{\linkS4class{packedMatrix}}, and
  \code{\linkS4class{sparseMatrix}},
  the specific class is
  \code{\linkS4class{Cholesky}},
  \code{\linkS4class{pCholesky}}, and
  \code{\linkS4class{dCHMsimpl}} or \code{\linkS4class{dCHMsuper}},
  respectively.
}
\details{
  Note that the result of a call to \code{Cholesky} inherits
  from \code{\linkS4class{CholeskyFactorization}} but not
  \code{\linkS4class{Matrix}}.  Users who just want a matrix
  should consider using \code{\link{chol}}, whose methods are
  simple wrappers around \code{Cholesky} returning just the
  upper triangular Cholesky factor \eqn{L'},
  typically as a \code{\linkS4class{triangularMatrix}}.
  However, a more principled approach would be to construct
  factors as needed from the \code{CholeskyFactorization} object,
  e.g., with \code{\link{expand1}(x, "L")}, if \code{x} is the
  object.

  The behaviour of \code{Cholesky(A, perm = TRUE)} for dense \code{A}
  is somewhat exceptional, in that it expects \emph{without} checking
  that \code{A} is positive semidefinite.  By construction, if \eqn{A}
  is positive semidefinite and the exact algorithm encounters a zero
  pivot, then the unfactorized trailing submatrix is the zero matrix,
  and there is nothing left to do.  Hence when the finite precision
  algorithm encounters a pivot less than \code{tol}, it signals a
  warning instead of an error and zeros the trailing submatrix in
  order to guarantee that \eqn{P' L L' P}{P' * L * L' * P} is positive
  semidefinite even if \eqn{A} is not.  It follows that one way to
  test for positive semidefiniteness of \eqn{A} in the event of a
  warning is to analyze the error
  \deqn{\frac{\lVert A - P' L L' P \rVert}{\lVert A \rVert}\,.}{
    norm(A - P' * L * L' * P) / norm(A).}
  See the examples and LAPACK Working Note (\dQuote{LAWN}) 161
  for details.
}
\seealso{
  Classes \code{\linkS4class{Cholesky}}, \code{\linkS4class{pCholesky}},
  \code{\linkS4class{dCHMsimpl}} and \code{\linkS4class{dCHMsuper}}
  and their methods.

  Classes \code{\linkS4class{dpoMatrix}}, \code{\linkS4class{dppMatrix}},
  and \code{\linkS4class{dsCMatrix}}.

  Generic function \code{\link{chol}},
  for obtaining the upper triangular Cholesky factor \eqn{L'} as a
  matrix or \code{\linkS4class{Matrix}}.

  Generic functions \code{\link{expand1}} and \code{\link{expand2}},
  for constructing matrix factors from the result.

  Generic functions \code{\link{BunchKaufman}}, \code{\link{Schur}},
  \code{\link{lu}}, and \code{\link{qr}},
  for computing other factorizations.
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

  Lucas, C. (2004).
  \emph{LAPACK-style codes for level 2 and 3 pivoted Cholesky factorizations}.
  LAPACK Working Note, Number 161.
  \url{https://www.netlib.org/lapack/lawnspdf/lawn161.pdf}
  
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
showMethods("Cholesky", inherited = FALSE)
set.seed(0)

## ---- Dense ----------------------------------------------------------

## .... Positive definite ..............................................

n <- 6L
(A1 <- crossprod(Matrix(rnorm(n * n), n, n)))
(ch.A1.nopivot <- Cholesky(A1, perm = FALSE))
(ch.A1 <- Cholesky(A1))
stopifnot(exprs = {
    length(ch.A1@perm) == ncol(A1)
    isPerm(ch.A1@perm)
    is.unsorted(ch.A1@perm) # typically not the identity permutation
    length(ch.A1.nopivot@perm) == 0L
})

## A ~ P1' L D L' P1 ~ P1' L L' P1 in floating point
str(e.ch.A1 <- expand2(ch.A1, LDL =  TRUE), max.level = 2L)
str(E.ch.A1 <- expand2(ch.A1, LDL = FALSE), max.level = 2L)
stopifnot(exprs = {
    all.equal(as(A1, "matrix"), as(Reduce(`\%*\%`, e.ch.A1), "matrix"))
    all.equal(as(A1, "matrix"), as(Reduce(`\%*\%`, E.ch.A1), "matrix"))
})

## .... Positive semidefinite but not positive definite ................

A2 <- A1
A2[1L, ] <- A2[, 1L] <- 0
A2
try(Cholesky(A2, perm = FALSE)) # fails as not positive definite
ch.A2 <- Cholesky(A2) # returns, with a warning and ...
A2.hat <- Reduce(`\%*\%`, expand2(ch.A2, LDL = FALSE))
norm(A2 - A2.hat, "2") / norm(A2, "2") # 7.670858e-17

## .... Not positive semidefinite ......................................

A3 <- A1
A3[1L, ] <- A3[, 1L] <- -1
A3
try(Cholesky(A3, perm = FALSE)) # fails as not positive definite
ch.A3 <- Cholesky(A3) # returns, with a warning and ...
A3.hat <- Reduce(`\%*\%`, expand2(ch.A3, LDL = FALSE))
norm(A3 - A3.hat, "2") / norm(A3, "2") # 1.781568

## Indeed, 'A3' is not positive semidefinite, but 'A3.hat' _is_
ch.A3.hat <- Cholesky(A3.hat)
A3.hat.hat <- Reduce(`\%*\%`, expand2(ch.A3.hat, LDL = FALSE))
norm(A3.hat - A3.hat.hat, "2") / norm(A3.hat, "2") # 1.777944e-16

## ---- Sparse ---------------------------------------------------------

## Really just three cases modulo permutation :
##
##            type        factorization  minors of P1 A P1'
##   1  simplicial  P1 A P1' = L1 D L1'             nonzero
##   2  simplicial  P1 A P1' = L    L '            positive
##   3  supernodal  P1 A P2' = L    L '            positive

data(KNex, package = "Matrix")
A4 <- crossprod(KNex[["mm"]])

ch.A4 <-
list(pivoted =
     list(simpl1 = Cholesky(A4, perm =  TRUE, super = FALSE, LDL =  TRUE),
          simpl0 = Cholesky(A4, perm =  TRUE, super = FALSE, LDL = FALSE),
          super0 = Cholesky(A4, perm =  TRUE, super =  TRUE             )),
     unpivoted =
     list(simpl1 = Cholesky(A4, perm = FALSE, super = FALSE, LDL =  TRUE),
          simpl0 = Cholesky(A4, perm = FALSE, super = FALSE, LDL = FALSE),
          super0 = Cholesky(A4, perm = FALSE, super =  TRUE             )))
ch.A4

s <- simplify2array
rapply2 <- function(object, f, ...) rapply(object, f, , , how = "list", ...)

s(rapply2(ch.A4, isLDL))
s(m.ch.A4 <- rapply2(ch.A4, expand1, "L")) # giving L = L1 sqrt(D)

## By design, the pivoted and simplicial factorizations
## are more sparse than the unpivoted and supernodal ones ...
s(rapply2(m.ch.A4, object.size))

## Which is nicely visualized by lattice-based methods for 'image'
inm <- c("pivoted", "unpivoted")
jnm <- c("simpl1", "simpl0", "super0")
for(i in 1:2)
for(j in 1:3)
print(image(m.ch.A4[[c(i, j)]], main = paste(inm[i], jnm[j])),
            split = c(j, i, 3L, 2L), more = i * j < 6L)

simpl1 <- ch.A4[[c("pivoted", "simpl1")]]
stopifnot(exprs = {
    length(simpl1@perm) == ncol(A4)
    isPerm(simpl1@perm, 0L)
    is.unsorted(simpl1@perm) # typically not the identity permutation
})

## One can expand with and without D regardless of isLDL(.),
## but "without" requires L = L1 sqrt(D), which is conditional
## on min(diag(D)) >= 0, hence "with" is the default
isLDL(simpl1)
stopifnot(min(diag(simpl1)) >= 0)
str(e.ch.A4 <- expand2(simpl1, LDL =  TRUE), max.level = 2L) # default
str(E.ch.A4 <- expand2(simpl1, LDL = FALSE), max.level = 2L)
stopifnot(exprs = {
    all.equal(E.ch.A4[["L" ]], e.ch.A4[["L1" ]] \%*\% sqrt(e.ch.A4[["D"]]))
    all.equal(E.ch.A4[["L."]], sqrt(e.ch.A4[["D"]]) \%*\% e.ch.A4[["L1."]])
    all.equal(A4, as(Reduce(`\%*\%`, e.ch.A4), "symmetricMatrix"))
    all.equal(A4, as(Reduce(`\%*\%`, E.ch.A4), "symmetricMatrix"))
})

## The "same" permutation matrix with "alternate" representation
## [i, perm[i]] {margin=1} <-> [invertPerm(perm)[j], j] {margin=2}
alt <- function(P) {
    P@margin <- 1L + !(P@margin - 1L) # 1 <-> 2
    P@perm <- invertPerm(P@perm)
    P
}

## Expansions are elegant but inefficient (transposes are redundant)
## hence programmers should consider methods for 'expand1' and 'diag'
stopifnot(exprs = {
    identical(expand1(simpl1, "P1"), alt(e.ch.A4[["P1"]]))
    identical(expand1(simpl1, "L"), E.ch.A4[["L"]])
    identical(Diagonal(x = diag(simpl1)), e.ch.A4[["D"]])
})

## chol(A, pivot = value) is a simple wrapper around
## Cholesky(A, perm = value, LDL = FALSE, super = FALSE),
## returning L' = sqrt(D) L1' _but_ giving no information
## about the permutation P1
selectMethod("chol", "dsCMatrix")
stopifnot(all.equal(chol(A4, pivot = TRUE), E.ch.A4[["L."]]))

## Now a symmetric matrix with positive _and_ negative eigenvalues,
## hence _not_ positive semidefinite
A5 <- new("dsCMatrix",
          Dim = c(7L, 7L),
          p = c(0:1, 3L, 6:7, 10:11, 15L),
          i = c(0L, 0:1, 0:3, 2:5, 3:6),
          x = c(1, 6, 38, 10, 60, 103, -4, 6, -32, -247, -2, -16, -128, -2, -67))
(ev <- eigen(A5, only.values = TRUE)$values)
(t.ev <- table(factor(sign(ev), -1:1))) # the matrix "inertia"

ch.A5 <- Cholesky(A5)
isLDL(ch.A5)
(d.A5 <- diag(ch.A5)) # diag(D) is partly negative

## Sylvester's law of inertia holds here, but not in general
## in finite precision arithmetic
stopifnot(identical(table(factor(sign(d.A5), -1:1)), t.ev))

try(expand1(ch.A5, "L"))         # unable to compute L = L1 sqrt(D)
try(expand2(ch.A5, LDL = FALSE)) # ditto
try(chol(A5, pivot = TRUE))      # ditto

## The default expansion is "square root free" and still works here
str(e.ch.A5 <- expand2(ch.A5, LDL = TRUE), max.level = 2L)
stopifnot(all.equal(A5, as(Reduce(`\%*\%`, e.ch.A5), "symmetricMatrix")))

## Version of the SuiteSparse library, which includes CHOLMOD
.SuiteSparse_version()
}
