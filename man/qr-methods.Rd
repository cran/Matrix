\name{qr-methods}
\title{Methods for QR Factorization}
%
\docType{methods}
\keyword{algebra}
\keyword{array}
\keyword{methods}
%
\alias{qr}
\alias{qr-methods}
%
\alias{qr,dgCMatrix-method}
\alias{qr,sparseMatrix-method}
%
\description{
  Computes the pivoted QR factorization of an \eqn{m \times n}{m-by-n}
  real matrix \eqn{A}, which has the general form
  \deqn{P_{1} A P_{2} = Q R}{P1 * A * P2 = Q * R}
  or (equivalently)
  \deqn{A = P_{1}' Q R P_{2}'}{A = P1' * Q * R * P2'}
  where
  \eqn{P_{1}}{P1} and \eqn{P_{2}}{P2} are permutation matrices,
  \eqn{Q = \prod_{j = 1}^{n} H_{j}}{Q = prod(Hj : j = 1,...,n)}
  is an \eqn{m \times m}{m-by-m} orthogonal matrix
  equal to the product of \eqn{n} Householder matrices \eqn{H_{j}}{Hj}, and
  \eqn{R} is an \eqn{m \times n}{m-by-n} upper trapezoidal matrix.
  
  \code{\linkS4class{denseMatrix}} use the default method implemented
  in \pkg{base}, namely \code{\link{qr.default}}.  It is built on
  LINPACK routine \code{dqrdc} and LAPACK routine \code{dgeqp3}, which
  do not pivot rows, so that \eqn{P_{1}}{P1} is an identity matrix.
  
  Methods for \code{\linkS4class{sparseMatrix}} are built on
  CSparse routines \code{cs_sqr} and \code{cs_qr}, which require
  \eqn{m \ge n}{m >= n}.
}
\usage{
qr(x, \dots)
\S4method{qr}{dgCMatrix}(x, order = 3L, \dots)
}
\arguments{
  \item{x}{a \link[=is.finite]{finite} matrix or
    \code{\linkS4class{Matrix}} to be factorized,
    satisfying \code{nrow(x) >= ncol(x)} if sparse.}
  \item{order}{an integer in \code{0:3} passed to CSparse routine
    \code{cs_sqr}, indicating a strategy for choosing the column
    permutation \eqn{P_{2}}{P2}.  0 means no column permutation.
    1, 2, and 3 indicate a fill-reducing ordering of \eqn{A + A'},
    \eqn{\tilde{A}' \tilde{A}}{A~' * A~}, and \eqn{A' A}{A' * A},
    where \eqn{\tilde{A}}{A~} is \eqn{A} with \dQuote{dense} rows
    removed.
    Do not set to 0 unless you know that the column order of \eqn{A}
    is already sensible.}
  \item{\dots}{further arguments passed to or from methods.}
}
\value{
  An object representing the factorization, inheriting from
  virtual S4 class \code{\linkS4class{QR}} or S3 class
  \code{\link[base]{qr}}.  The specific class is \code{qr}
  unless \code{x} inherits from virtual class
  \code{\linkS4class{sparseMatrix}}, in which case it is
  \code{\linkS4class{sparseQR}}.
}
\details{
  If \code{x} is sparse and structurally rank deficient, having
  structural rank \eqn{r < n}, then \code{x} is augmented with
  \eqn{(n-r)} rows of (partly non-structural) zeros, such that
  the augmented matrix has structural rank \eqn{n}.
  This augmented matrix is factorized as described above:
  \Sdeqn{P_1 A P_2 = P_1 \\\\\\\\begin{bmatrix} A_{0} \\\\\\\\\\\\\\\\ 0 \\\\\\\\end{bmatrix} P_2 = Q R}{P1 * A * P2 = P1 * [A0; 0] * P2 = Q * R}
  where \eqn{A_0}{A0} denotes the original, user-supplied
  \eqn{(m-(n-r)) \times n}{(m-(n-r))-by-n} matrix.
}
\seealso{
  Class \code{\linkS4class{sparseQR}} and its methods.

  Class \code{\linkS4class{dgCMatrix}}.
  
  Generic function \code{\link[base]{qr}} from \pkg{base},
  whose default method \code{qr.default} \dQuote{defines}
  the S3 class \code{qr} of dense QR factorizations.

  Generic functions \code{\link{expand1}} and \code{\link{expand2}},
  for constructing matrix factors from the result.
  
  Generic functions \code{\link{Cholesky}}, \code{\link{BunchKaufman}},
  \code{\link{Schur}}, and \code{\link{lu}},
  for computing other factorizations.
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
showMethods("qr", inherited = FALSE)

## Rank deficient: columns 3 {b2} and 6 {c3} are "extra"
M <- as(cbind(a1 = 1,
              b1 = rep(c(1, 0), each = 3L),
              b2 = rep(c(0, 1), each = 3L),
              c1 = rep(c(1, 0, 0), 2L),
              c2 = rep(c(0, 1, 0), 2L),
              c3 = rep(c(0, 0, 1), 2L)),
        "CsparseMatrix")
rownames(M) <- paste0("r", seq_len(nrow(M)))
b <- 1:6
eps <- .Machine$double.eps

## .... [1] full rank ..................................................
## ===> a least squares solution of A x = b exists
##      and is unique _in exact arithmetic_

(A1 <- M[, -c(3L, 6L)])
(qr.A1 <- qr(A1))

stopifnot(exprs = {
    rankMatrix(A1) == ncol(A1)
    { d1 <- abs(diag(qr.A1@R)); sum(d1 < max(d1) * eps) == 0L }
    rcond(crossprod(A1)) >= eps
    all.equal(qr.coef(qr.A1, b), drop(solve(crossprod(A1), crossprod(A1, b))))
    all.equal(qr.fitted(qr.A1, b) + qr.resid(qr.A1, b), b)
})

## .... [2] numerically rank deficient with full structural rank .......
## ===> a least squares solution of A x = b does not
##      exist or is not unique _in exact arithmetic_

(A2 <- M)
(qr.A2 <- qr(A2))

stopifnot(exprs = {
    rankMatrix(A2) == ncol(A2) - 2L
    { d2 <- abs(diag(qr.A2@R)); sum(d2 < max(d2) * eps) == 2L }
    rcond(crossprod(A2)) < eps

    ## 'qr.coef' computes unique least squares solution of "nearby" problem
    ## Z x = b for some full rank Z ~ A, currently without warning {FIXME} !
    tryCatch({ qr.coef(qr.A2, b); TRUE }, condition = function(x) FALSE)

    all.equal(qr.fitted(qr.A2, b) + qr.resid(qr.A2, b), b)
})

## .... [3] numerically and structurally rank deficient ................
## ===> factorization of _augmented_ matrix with
##      full structural rank proceeds as in [2]

##  NB: implementation details are subject to change; see (*) below

A3 <- M
A3[, c(3L, 6L)] <- 0
A3
(qr.A3 <- qr(A3)) # with a warning ... "additional 2 row(s) of zeros"

stopifnot(exprs = {
    ## sparseQR object preserves the unaugmented dimensions (*)
    dim(qr.A3  ) == dim(A3)
    dim(qr.A3@V) == dim(A3) + c(2L, 0L)
    dim(qr.A3@R) == dim(A3) + c(2L, 0L)

    ## The augmented matrix remains numerically rank deficient
    rankMatrix(A3) == ncol(A3) - 2L
    { d3 <- abs(diag(qr.A3@R)); sum(d3 < max(d3) * eps) == 2L }
    rcond(crossprod(A3)) < eps
})

## Auxiliary functions accept and return a vector or matrix
## with dimensions corresponding to the unaugmented matrix (*),
## in all cases with a warning
qr.coef  (qr.A3, b)
qr.fitted(qr.A3, b)
qr.resid (qr.A3, b)

## .... [4] yet more examples ..........................................

## By disabling column pivoting, one gets the "vanilla" factorization
## A = Q~ R, where Q~ := P1' Q is orthogonal because P1 and Q are

(qr.A1.pp <- qr(A1, order = 0L)) # partial pivoting

ae1 <- function(a, b, ...) all.equal(as(a, "matrix"), as(b, "matrix"), ...)
ae2 <- function(a, b, ...) ae1(unname(a), unname(b), ...)

stopifnot(exprs = {
    length(qr.A1   @q) == ncol(A1)
    length(qr.A1.pp@q) == 0L # indicating no column pivoting
    ae2(A1[, qr.A1@q + 1L], qr.Q(qr.A1   ) \%*\% qr.R(qr.A1   ))
    ae2(A1                , qr.Q(qr.A1.pp) \%*\% qr.R(qr.A1.pp))
})
}
