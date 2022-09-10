\name{solve-methods}
\title{Methods in Package Matrix for Function \code{solve()}}
\docType{methods}
%
\alias{solve} % catch ?solve too
\alias{solve-methods}
%
\alias{solve,ANY,ANY-method}
\alias{solve,CHMfactor,denseMatrix-method}
\alias{solve,CHMfactor,sparseMatrix-method}
\alias{solve,CHMfactor,matrix-method}
\alias{solve,CHMfactor,missing-method}
\alias{solve,CHMfactor,numLike-method}
\alias{solve,CsparseMatrix,ANY-method}
\alias{solve,Matrix,sparseVector-method}
\alias{solve,MatrixFactorization,ANY-method}
\alias{solve,MatrixFactorization,missing-method}
\alias{solve,MatrixFactorization,sparseVector-method}
\alias{solve,RsparseMatrix,ANY-method}
\alias{solve,TsparseMatrix,ANY-method}
\alias{solve,ddiMatrix,Matrix-method}
\alias{solve,ddiMatrix,matrix-method}
\alias{solve,ddiMatrix,missing-method}
\alias{solve,ddiMatrix,numLike-method}
\alias{solve,denseLU,missing-method}
\alias{solve,denseMatrix,ANY-method}
\alias{solve,dgCMatrix,denseMatrix-method}
\alias{solve,dgCMatrix,matrix-method}
\alias{solve,dgCMatrix,missing-method}
\alias{solve,dgCMatrix,numLike-method}
\alias{solve,dgCMatrix,sparseMatrix-method}
\alias{solve,dgeMatrix,Matrix-method}
\alias{solve,dgeMatrix,matrix-method}
\alias{solve,dgeMatrix,missing-method}
\alias{solve,dgeMatrix,numLike-method}
\alias{solve,diagonalMatrix,ANY-method}
\alias{solve,dpoMatrix,Matrix-method}
\alias{solve,dpoMatrix,matrix-method}
\alias{solve,dpoMatrix,missing-method}
\alias{solve,dpoMatrix,numLike-method}
\alias{solve,dppMatrix,Matrix-method}
\alias{solve,dppMatrix,matrix-method}
\alias{solve,dppMatrix,missing-method}
\alias{solve,dppMatrix,numLike-method}
\alias{solve,dsCMatrix,denseMatrix-method}
\alias{solve,dsCMatrix,matrix-method}
\alias{solve,dsCMatrix,missing-method}
\alias{solve,dsCMatrix,numLike-method}
\alias{solve,dsCMatrix,sparseMatrix-method}
\alias{solve,dspMatrix,Matrix-method}
\alias{solve,dspMatrix,matrix-method}
\alias{solve,dspMatrix,missing-method}
\alias{solve,dspMatrix,numLike-method}
\alias{solve,dsyMatrix,Matrix-method}
\alias{solve,dsyMatrix,matrix-method}
\alias{solve,dsyMatrix,missing-method}
\alias{solve,dsyMatrix,numLike-method}
\alias{solve,dtpMatrix,Matrix-method}
\alias{solve,dtpMatrix,matrix-method}
\alias{solve,dtpMatrix,missing-method}
\alias{solve,dtpMatrix,numLike-method}
\alias{solve,dtrMatrix,Matrix-method}
\alias{solve,dtrMatrix,matrix-method}
\alias{solve,dtrMatrix,missing-method}
\alias{solve,dtrMatrix,numLike-method}
\alias{solve,indMatrix,ANY-method}
\alias{solve,pMatrix,Matrix-method}
\alias{solve,pMatrix,matrix-method}
\alias{solve,pMatrix,missing-method}
\alias{solve,pMatrix,numLike-method}
\alias{solve,dtCMatrix,dgCMatrix-method}
\alias{solve,dtCMatrix,dgeMatrix-method}
\alias{solve,dtCMatrix,dsCMatrix-method}
\alias{solve,dtCMatrix,dspMatrix-method}
\alias{solve,dtCMatrix,dsyMatrix-method}
\alias{solve,dtCMatrix,dtCMatrix-method}
\alias{solve,dtCMatrix,dtpMatrix-method}
\alias{solve,dtCMatrix,dtrMatrix-method}
\alias{solve,dtCMatrix,denseMatrix-method}
\alias{solve,dtCMatrix,matrix-method}
\alias{solve,dtCMatrix,missing-method}
\alias{solve,dtCMatrix,numLike-method}
\alias{solve,dtCMatrix,sparseMatrix-method}
\alias{solve,matrix,Matrix-method}
\alias{solve,matrix,sparseVector-method}
\alias{solve,sparseQR,ANY-method}
\alias{solve,sparseQR,missing-method}
%
\description{
  Methods for generic function \code{\link{solve}}, for solving linear
  systems of equations.  These solve for \eqn{X} in \deqn{A X = B}
  where \eqn{A} is a square matrix and \eqn{X} and \eqn{B} are matrices
  with compatible dimensions.  The usual \R syntax is
  \preformatted{
                        x <- solve(a, b, ...)
  }
  where \code{b} may also be a vector, in which case it is treated
  as a 1-column matrix.  Methods support \code{a} inheriting from
  virtual classes
  \code{\linkS4class{Matrix}} and \code{\linkS4class{MatrixFactorization}}
  and \code{b} inheriting from virtual classes
  \code{Matrix} and \code{\linkS4class{sparseVector}}.
}
\usage{
## solve(a, b, ...) # the two-argument version, almost always preferred to
## solve(a,    ...) # the *rarely needed* one-argument version
%^must comment; otherwise 'R CMD check' waffles
\S4method{solve}{dgCMatrix,missing}(a, b, sparse = NA, \dots)
\S4method{solve}{dgCMatrix,matrix}(a, b, sparse = FALSE, \dots)
\S4method{solve}{dgCMatrix,denseMatrix}(a, b, sparse = FALSE, \dots)
\S4method{solve}{dgCMatrix,sparseMatrix}(a, b, sparse = NA, tol = .Machine$double.eps, \dots)
\S4method{solve}{CHMfactor,denseMatrix}(a, b, system = c("A", "LDLt", "LD", "DLt", "L", "Lt", "D", "P", "Pt"), \dots)
}
\arguments{
  \item{a}{a square numeric matrix, \eqn{A}, typically of one of the
    classes in \pkg{Matrix}.  Logical matrices are coerced to
    corresponding numeric ones.}
  \item{b}{numeric vector or matrix (dense or sparse) as RHS
    of the linear system \eqn{Ax = b}.}
  \item{sparse}{only when \code{a} is a \code{\linkS4class{sparseMatrix}}: 
    logical specifying if the result should also (formally) be sparse.} % FIXME: mention pos. def., etc.
  \item{tol}{only when \code{a} is a \code{sparseMatrix}
    and \code{sparse} is \code{TRUE}: an error is signaled
    if the ratio \code{min(d)/max(d)}, where \code{d = abs(diag(U))}
    and \code{A = L U}, is \emph{less} than \code{tol}, indicating
    (near-)singular \eqn{A}.}
%  \item{tol}{only used when \code{a} is sparse, in the
%      \code{\link{isSymmetric}(a, tol=*)} test, where that applies.}
  \item{system}{only when \code{a} is a \code{\linkS4class{CHMfactor}}:
    character string indicating the kind of linear system to be
    solved, see below.  Note that the default, \code{"A"}, does
    \emph{not} solve the triangular system (but \code{"L"} does).}
  \item{\dots}{potentially further arguments to the methods.}
}

\section{Methods}{
  \describe{
    \item{\code{signature(a = "ANY", b = "ANY")}}{is simply the
      \pkg{base} package's S3 generic \code{\link{solve}}.}

%% This is copy-paste  in  CHMfactor-class.Rd  {FIXME ?}
    \item{\code{signature(a = "CHMfactor", b = "...."), system= *}}{The
      \code{solve} methods for a \code{"\linkS4class{CHMfactor}"} object
      take an optional third argument \code{system} whose value can be
      one of the character strings \code{"A"}, \code{"LDLt"}, \code{"LD"},
      \code{"DLt"}, \code{"L"}, \code{"Lt"}, \code{"D"}, \code{"P"} or
      \code{"Pt"}.  This argument describes the system to be solved.  The
      default,  \code{"A"}, is to solve \eqn{Ax = b} for \eqn{x} where
      \code{A} is sparse, positive-definite matrix that was factored to produce
      \code{a}.  Analogously, \code{system = "L"} returns the solution
      \eqn{x}, of \eqn{Lx = b}; similarly, for all system codes
      \bold{but} \code{"P"} and \code{"Pt"} where, e.g., \code{x <-
	solve(a, b,system="P")} is equivalent to \code{x <- P \%*\% b}.

      If \code{b} is a \code{\linkS4class{sparseMatrix}}, \code{system}
      is used as above the corresponding sparse CHOLMOD algorithm is called.
    }

    \item{\code{signature(a = "ddenseMatrix", b = "....")}}{(for all
      \code{b}) work via \code{as(a, "generalMatrix")}, using the its
      methods, see below.}

    \item{\code{signature(a = "denseLU", b = "missing")}}{
      basically computes uses triangular forward- and back-solve.}

    \item{\code{signature(a = "dgCMatrix", b = "matrix")}}{, and}
    %% -> ../R/dgCMatrix.R
    \item{\code{signature(a = "dgCMatrix", b = "ddenseMatrix")}}{with extra
      argument list \code{( sparse = FALSE, tol = .Machine$double.eps ) }:
      Uses the sparse \code{\link{lu}(a)} decomposition (which is cached
      in \code{a}'s \code{factor} slot).
      By default, \code{sparse=FALSE}, returns a
      \code{\linkS4class{denseMatrix}}, since \eqn{U^{-1} L^{-1} B} may
      not be sparse at all, even when \eqn{L} and \eqn{U} are.

      If \code{sparse=TRUE}, returns a \code{\linkS4class{sparseMatrix}}
      (which may not be very sparse at all, even if \code{a} \emph{was} sparse).
    }

    \item{\code{signature(a = "dgCMatrix", b = "dsparseMatrix")}}{, and}
    \item{\code{signature(a = "dgCMatrix", b = "missing")}}{with extra
      argument list \code{( sparse=FALSE, tol = .Machine$double.eps ) }:
      Checks if \code{a} is symmetric, and in that case, coerces it to
      \code{"\linkS4class{symmetricMatrix}"}, and then computes a
      \emph{sparse} solution via sparse Cholesky factorization,
      independently of the \code{sparse} argument.  If \code{a} is not
      symmetric, the sparse \code{\link{lu}} decomposition is used
      and the result will be sparse or dense, depending on the
      \code{sparse} argument, exactly as for the above (\code{b =
	"ddenseMatrix"}) case.
    }

    \item{\code{signature(a = "dgeMatrix", b = ".....")}}{
      solve the system via internal LU, calling LAPACK routines
      \code{dgetri} or \code{dgetrs}.
    }

    \item{\code{signature(a = "diagonalMatrix", b = "matrix")}}{and
      other \code{b}s: Of course this is trivially implemented, as
      \eqn{D^{-1}} is diagonal with entries \eqn{1 / D[i,i]}.}

    \item{\code{signature(a = "dpoMatrix", b = "....Matrix")}}{, and}
    \item{\code{signature(a = "dppMatrix", b = "....Matrix")}}{
      The Cholesky decomposition of \code{a} is calculated (if
      needed) while solving the system.}

    \item{\code{signature(a = "dsCMatrix", b = "....")}}{% ../R/dsCMatrix.R
      All these methods first try Cholmod's Cholesky factorization; if
      that works, i.e., typically if \code{a} is positive semi-definite,
      it is made use of.  Otherwise, the sparse LU decomposition is used
      as for the \dQuote{general} matrices of class \code{"dgCMatrix"}.}

    \item{\code{signature(a = "dspMatrix", b = "....")}}{, and}
    \item{\code{signature(a = "dsyMatrix", b = "....")}}{% ../R/dsyMatrix.R
      all end up calling LAPACK routines \code{dsptri}, \code{dsptrs},
      \code{dsytrs} and \code{dsytri}.
    }

    \item{\code{signature(a = "dtCMatrix", b = "CsparseMatrix")}}{,}
    \item{\code{signature(a = "dtCMatrix", b = "dgeMatrix")}}{, etc
      sparse triangular solve, in traditional S/\R also known as
      \code{\link{backsolve}}, or \code{\link{forwardsolve}}.
      \code{solve(a,b)} is a \code{\linkS4class{sparseMatrix}} if
      \code{b} is, and hence a \code{\linkS4class{denseMatrix}}
      otherwise.
    }

    \item{\code{signature(a = "dtrMatrix", b = "ddenseMatrix")}}{, and}
    \item{\code{signature(a = "dtpMatrix", b = "matrix")}}{, and% ../R/dtrMatrix.R
      similar \code{b}, including \code{"missing"}, and
      \code{"diagonalMatrix"}:

      all use LAPACK based versions of efficient triangular
      \code{\link{backsolve}}, or \code{\link{forwardsolve}}.
    }

    \item{\code{signature(a = "Matrix", b = "diagonalMatrix")}}{
      works via \code{as(b, "CsparseMatrix")}.}

    \item{\code{signature(a = "sparseQR", b = "ANY")}}{
      simply uses \code{\link{qr.coef}(a, b)}.}

    \item{\code{signature(a = "pMatrix", b = ".....")}}{
      these methods typically use \code{\link{crossprod}(a,b)}, as
      the inverse of a permutation matrix is the same as its transpose.}

    \item{\code{signature(a = "TsparseMatrix", b = "ANY")}}{
      all work via \code{as(a, "CsparseMatrix")}.}

  }
}%{Methods}
\seealso{
  \code{\link{solve}}, \code{\link{lu}}, and class documentations
  \code{\linkS4class{CHMfactor}}, \code{\linkS4class{sparseLU}}, and
  \code{\linkS4class{MatrixFactorization}}.
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
image(ia0 <- solve(a)) # checker board, dense [but really, a is singular!]
try(solve(a, sparse=TRUE))##-> error [ TODO: assertError ]
ia. <- solve(a, sparse=TRUE, tol = 1e-19)##-> *no* error
if(R.version$arch == "x86_64")
  ## Fails on 32-bit [Fedora 19, R 3.0.2] from Matrix 1.1-0 on [FIXME ??] only
  stopifnot(all.equal(as.matrix(ia.), as.matrix(ia0)))
a <- a + Diagonal(n)
iad <- solve(a)
ias <- solve(a, sparse=TRUE)
stopifnot(all.equal(as(ias,"denseMatrix"), iad, tolerance=1e-14))
I. <- iad \%*\% a          ; image(I.)
I0 <- drop0(zapsmall(I.)); image(I0)
.I <- a \%*\% iad
.I0 <- drop0(zapsmall(.I))
stopifnot( all.equal(as(I0, "diagonalMatrix"), Diagonal(n)),
           all.equal(as(.I0,"diagonalMatrix"), Diagonal(n)) )

}
\keyword{methods}
