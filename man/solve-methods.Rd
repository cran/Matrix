\name{solve-methods}
\title{Methods in Package \pkg{Matrix} for Function \code{solve}}
%
\docType{methods}
\keyword{algebra}
\keyword{array}
\keyword{methods}
%
\alias{solve}
\alias{solve-methods}
%
\alias{solve,ANY,ANY-method}
\alias{solve,BunchKaufman,missing-method}
\alias{solve,BunchKaufman,dgeMatrix-method}
\alias{solve,CHMfactor,missing-method}
\alias{solve,CHMfactor,dgeMatrix-method}
\alias{solve,CHMfactor,dgCMatrix-method}
\alias{solve,Cholesky,missing-method}
\alias{solve,Cholesky,dgeMatrix-method}
\alias{solve,CsparseMatrix,ANY-method}
\alias{solve,Matrix,sparseVector-method}
\alias{solve,MatrixFactorization,CsparseMatrix-method}
\alias{solve,MatrixFactorization,RsparseMatrix-method}
\alias{solve,MatrixFactorization,TsparseMatrix-method}
\alias{solve,MatrixFactorization,denseMatrix-method}
\alias{solve,MatrixFactorization,dgCMatrix-method}
\alias{solve,MatrixFactorization,dgeMatrix-method}
\alias{solve,MatrixFactorization,diagonalMatrix-method}
\alias{solve,MatrixFactorization,indMatrix-method}
\alias{solve,MatrixFactorization,matrix-method}
\alias{solve,MatrixFactorization,sparseVector-method}
\alias{solve,MatrixFactorization,vector-method}
\alias{solve,RsparseMatrix,ANY-method}
\alias{solve,Schur,ANY-method}
\alias{solve,TsparseMatrix,ANY-method}
\alias{solve,ddiMatrix,Matrix-method}
\alias{solve,ddiMatrix,matrix-method}
\alias{solve,ddiMatrix,missing-method}
\alias{solve,ddiMatrix,vector-method}
\alias{solve,denseLU,missing-method}
\alias{solve,denseLU,dgeMatrix-method}
\alias{solve,denseMatrix,ANY-method}
\alias{solve,dgCMatrix,denseMatrix-method}
\alias{solve,dgCMatrix,matrix-method}
\alias{solve,dgCMatrix,missing-method}
\alias{solve,dgCMatrix,sparseMatrix-method}
\alias{solve,dgCMatrix,vector-method}
\alias{solve,dgeMatrix,ANY-method}
\alias{solve,diagonalMatrix,ANY-method}
\alias{solve,dpoMatrix,ANY-method}
\alias{solve,dppMatrix,ANY-method}
\alias{solve,dsCMatrix,denseMatrix-method}
\alias{solve,dsCMatrix,matrix-method}
\alias{solve,dsCMatrix,missing-method}
\alias{solve,dsCMatrix,sparseMatrix-method}
\alias{solve,dsCMatrix,vector-method}
\alias{solve,dspMatrix,ANY-method}
\alias{solve,dsyMatrix,ANY-method}
\alias{solve,dtCMatrix,dgCMatrix-method}
\alias{solve,dtCMatrix,dgeMatrix-method}
\alias{solve,dtCMatrix,missing-method}
\alias{solve,dtCMatrix,triangularMatrix-method}
\alias{solve,dtpMatrix,dgeMatrix-method}
\alias{solve,dtpMatrix,missing-method}
\alias{solve,dtpMatrix,triangularMatrix-method}
\alias{solve,dtrMatrix,dgeMatrix-method}
\alias{solve,dtrMatrix,missing-method}
\alias{solve,dtrMatrix,triangularMatrix-method}
\alias{solve,indMatrix,ANY-method}
\alias{solve,matrix,Matrix-method}
\alias{solve,matrix,sparseVector-method}
\alias{solve,pBunchKaufman,missing-method}
\alias{solve,pBunchKaufman,dgeMatrix-method}
\alias{solve,pCholesky,missing-method}
\alias{solve,pCholesky,dgeMatrix-method}
\alias{solve,pMatrix,Matrix-method}
\alias{solve,pMatrix,matrix-method}
\alias{solve,pMatrix,missing-method}
\alias{solve,pMatrix,vector-method}
\alias{solve,sparseLU,missing-method}
\alias{solve,sparseLU,dgeMatrix-method}
\alias{solve,sparseLU,dgCMatrix-method}
\alias{solve,sparseQR,missing-method}
\alias{solve,sparseQR,dgeMatrix-method}
\alias{solve,sparseQR,dgCMatrix-method}
\alias{solve,triangularMatrix,CsparseMatrix-method}
\alias{solve,triangularMatrix,RsparseMatrix-method}
\alias{solve,triangularMatrix,TsparseMatrix-method}
\alias{solve,triangularMatrix,denseMatrix-method}
\alias{solve,triangularMatrix,dgCMatrix-method}
\alias{solve,triangularMatrix,dgeMatrix-method}
\alias{solve,triangularMatrix,diagonalMatrix-method}
\alias{solve,triangularMatrix,indMatrix-method}
\alias{solve,triangularMatrix,matrix-method}
\alias{solve,triangularMatrix,vector-method}
%
\description{
  Methods for generic function \code{\link[base]{solve}} for solving
  linear systems of equations,
  i.e., for \eqn{X} in \eqn{A X = B}{A * X = B},
  where \eqn{A} is a square matrix and \eqn{X} and \eqn{B} are matrices
  with dimensions consistent with \eqn{A}.
}
\usage{
solve(a, b, ...)

\S4method{solve}{dgeMatrix,ANY}(a, b, tol = .Machine$double.eps, \dots)

\S4method{solve}{dgCMatrix,missing}(a, b, sparse = TRUE, \dots)
\S4method{solve}{dgCMatrix,matrix}(a, b, sparse = FALSE, \dots)
\S4method{solve}{dgCMatrix,denseMatrix}(a, b, sparse = FALSE, \dots)
\S4method{solve}{dgCMatrix,sparseMatrix}(a, b, sparse = TRUE, \dots)

\S4method{solve}{denseLU,dgeMatrix}(a, b, \dots)
\S4method{solve}{BunchKaufman,dgeMatrix}(a, b, \dots)
\S4method{solve}{Cholesky,dgeMatrix}(a, b, \dots)
\S4method{solve}{sparseLU,dgCMatrix}(a, b, tol = .Machine$double.eps, \dots)
\S4method{solve}{sparseQR,dgCMatrix}(a, b, \dots)
\S4method{solve}{CHMfactor,dgCMatrix}(a, b, system = c("A", "LDLt", "LD", "DLt", "L", "Lt", "D", "P", "Pt"), \dots)
}
\arguments{
  \item{a}{a \link[=is.finite]{finite} square matrix or
    \code{\linkS4class{Matrix}} containing the coefficients
    of the linear system, or otherwise a
    \code{\linkS4class{MatrixFactorization}},
    in which case methods behave (by default)
    as if the factorized matrix were specified.}
  \item{b}{a vector, \code{\linkS4class{sparseVector}},
    matrix, or \code{\linkS4class{Matrix}} satisfying
    \code{NROW(b) == nrow(a)}, giving the right-hand side(s)
    of the linear system.  Vectors \code{b} are treated as
    \code{length(b)}-by-1 matrices.  If \code{b} is missing,
    then methods take \code{b} to be an identity matrix.}
  \item{tol}{a non-negative number.  For \code{a} inheriting from
    \code{\linkS4class{denseMatrix}}, an error is signaled if the
    reciprocal one-norm condition number (see \code{\link[base]{rcond}})
    of \code{a} is less than \code{tol}, indicating that \code{a} is
    near-singular.  For \code{a} of class \code{\linkS4class{sparseLU}},
    an error is signaled if the ratio \code{min(d)/max(d)} is less
    than \code{tol}, where \code{d = abs(diag(a@U))}.  (Interpret
    with care, as this ratio is a cheap heuristic and \emph{not}
    in general equal to or even proportional to the reciprocal
    one-norm condition number.)  Setting \code{tol = 0} disables
    the test.}
  \item{sparse}{a logical indicating if the result should be formally
    sparse, i.e., if the result should inherit from virtual class
    \code{\linkS4class{sparseMatrix}}.
    Only methods for sparse \code{a} and missing or matrix \code{b}
    have this argument.
    Methods for missing or sparse \code{b} use \code{sparse = TRUE}
    by default.  Methods for dense \code{b} use \code{sparse = FALSE}
    by default.}
  \item{system}{a string specifying a linear system to be solved.
    Only methods for \code{a}
    inheriting from \code{\linkS4class{CHMfactor}} have this argument.
    See \sQuote{Details}.}
  \item{\dots}{further arguments passed to or from methods.}
}
\details{
  Methods for general and symmetric matrices \code{a} compute a
  triangular factorization (LU, Bunch-Kaufman, or Cholesky)
  and call the method for the corresponding factorization class.
  The factorization is sparse if \code{a} is.  Methods for sparse,
  symmetric matrices \code{a} attempt a Cholesky factorization
  and perform an LU factorization only if that fails (typically
  because \code{a} is not positive definite).

  Triangular, diagonal, and permutation matrices do not require
  factorization (they are already \dQuote{factors}), hence methods
  for those are implemented directly.  For triangular \code{a},
  solutions are obtained by forward or backward substitution;
  for diagonal \code{a}, they are obtained by scaling the rows
  of \code{b}; and for permutations \code{a}, they are obtained
  by permuting the rows of \code{b}.

  Methods for dense \code{a} are built on 14 LAPACK routines:
  class \code{d..Matrix}, where \code{..=(ge|tr|tp|sy|sp|po|pp)},
  uses routines \code{d..tri} and \code{d..trs} for missing
  and non-missing \code{b}, respectively.  A corollary is that
  these methods always give a dense result.
  
  Methods for sparse \code{a} are built on CXSparse routines
  \code{cs_lsolve}, \code{cs_usolve}, and \code{cs_spsolve} and
  CHOLMOD routines \code{cholmod_solve} and \code{cholmod_spsolve}.
  By default, these methods give a vector result if \code{b}
  is a vector, a sparse matrix result if \code{b} is missing
  or a sparse matrix, and a dense matrix result if \code{b}
  is a dense matrix.  One can override this behaviour by setting
  the \code{sparse} argument, where available, but that should
  be done with care.  Note that a sparse result may be sparse only
  in the formal sense and not at all in the mathematical sense,
  depending on the nonzero patterns of \code{a} and \code{b}.
  Furthermore, whereas dense results are fully preallocated,
  sparse results must be \dQuote{grown} in a loop over the columns
  of \code{b}.
  
  Methods for \code{a} of class \code{\linkS4class{sparseQR}}
  are simple wrappers around \code{\link{qr.coef}}, giving the
  \emph{least squares} solution in overdetermined cases.
  
  Methods for \code{a} inheriting from \code{\linkS4class{CHMfactor}}
  can solve systems other than the default one \eqn{A X = B}{A * X = B}.
  The correspondence between its \code{system} argument the system
  actually solved is outlined in the table below.
  See \code{\link{CHMfactor-class}} for a definition of notation.
  
  \tabular{rrr}{
    \code{system} \tab \code{\link{isLDL}(a)=TRUE} \tab \code{\link{isLDL}(a)=FALSE}\cr
    \code{"A"} \tab \eqn{A X = B}{A * X = B} \tab \eqn{A X = B}{A * X = B}\cr
    \code{"LDLt"} \tab \eqn{L_{1} D L_{1}' X = B}{L1 * D * L1' * X = B} \tab \eqn{L L' X = B}{L * L' * X = B}\cr
    \code{"LD"} \tab \eqn{L_{1} D X = B}{L1 * D * X = B} \tab \eqn{L X = B}{L * X = B}\cr
    \code{"DLt"} \tab \eqn{D L_{1}' X = B}{D * L1' * X = B} \tab \eqn{L' X = B}{L' * X = B}\cr
    \code{"L"} \tab \eqn{L_{1} X = B}{L1 * X = B} \tab \eqn{L X = B}{L * X = B}\cr
    \code{"Lt"} \tab \eqn{L_{1}' X = B}{L1' * X = B} \tab \eqn{L' X = B}{L' * X = B}\cr
    \code{"D"} \tab \eqn{D X = B}{D * X = B} \tab \eqn{X = B}{X = B}\cr
    \code{"P"} \tab \eqn{X = P_{1} B}{X = P1 * B} \tab \eqn{X = P_{1} B}{X = P1 * B}\cr
    \code{"Pt"} \tab \eqn{X = P_{1}' B}{X = P1' * B} \tab \eqn{X = P_{1}' B}{X = P1' * B}
  }
}
\seealso{
  Virtual class \code{\linkS4class{MatrixFactorization}} and its
  subclasses.

  Generic functions \code{\link{Cholesky}}, \code{\link{BunchKaufman}},
  \code{\link{Schur}}, \code{\link{lu}}, and \code{\link{qr}} for
  \emph{computing} factorizations.

  Generic function \code{\link[base]{solve}} from \pkg{base}.

  Function \code{\link{qr.coef}} from \pkg{base} for computing
  least squares solutions of overdetermined linear systems.
}
\examples{
## A close to symmetric example with "quite sparse" inverse:
n1 <- 7; n2 <- 3
dd <- data.frame(a = gl(n1,n2), b = gl(n2,1,n1*n2))# balanced 2-way
X <- sparse.model.matrix(~ -1+ a + b, dd)# no intercept --> even sparser
XXt <- tcrossprod(X)
diag(XXt) <- rep(c(0,0,1,0), length.out = nrow(XXt))

n <- nrow(ZZ <- kronecker(XXt, Diagonal(x=c(4,1))))
image(a <- 2*Diagonal(n) + ZZ \%*\% Diagonal(x=c(10, rep(1, n-1))))
isSymmetric(a) # FALSE
image(drop0(skewpart(a)))
image(ia0 <- solve(a, tol = 0)) # checker board, dense [but really, a is singular!]
try(solve(a, sparse=TRUE))##-> error [ TODO: assertError ]
ia. <- solve(a, sparse=TRUE, tol = 1e-19)##-> *no* error
if(R.version$arch == "x86_64")
  ## Fails on 32-bit [Fedora 19, R 3.0.2] from Matrix 1.1-0 on [FIXME ??] only
  stopifnot(all.equal(as.matrix(ia.), as.matrix(ia0)))
a <- a + Diagonal(n)
iad <- solve(a)
ias <- solve(a, sparse=FALSE)
stopifnot(all.equal(as(iad,"denseMatrix"), ias, tolerance=1e-14))
I. <- iad \%*\% a          ; image(I.)
I0 <- drop0(zapsmall(I.)); image(I0)
.I <- a \%*\% iad
.I0 <- drop0(zapsmall(.I))
stopifnot( all.equal(as(I0, "diagonalMatrix"), Diagonal(n)),
           all.equal(as(.I0,"diagonalMatrix"), Diagonal(n)) )

}
