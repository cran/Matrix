\name{matrix-products}
\docType{methods}
\title{Matrix (Cross) Products (of Transpose)}
\alias{\%*\%-methods}
\alias{crossprod-methods}
\alias{tcrossprod-methods}
%%-- %*% -------------------------------------
\alias{\%*\%,ddenseMatrix,ldenseMatrix-method}
\alias{\%*\%,ddenseMatrix,matrix-method}
\alias{\%*\%,ddenseMatrix,ndenseMatrix-method}
\alias{\%*\%,ldenseMatrix,ddenseMatrix-method}
\alias{\%*\%,ldenseMatrix,ldenseMatrix-method}
\alias{\%*\%,ldenseMatrix,matrix-method}
\alias{\%*\%,ldenseMatrix,ndenseMatrix-method}
\alias{\%*\%,matrix,ddenseMatrix-method}
\alias{\%*\%,matrix,ldenseMatrix-method}
\alias{\%*\%,matrix,ndenseMatrix-method}
\alias{\%*\%,ndenseMatrix,ddenseMatrix-method}
\alias{\%*\%,ndenseMatrix,ldenseMatrix-method}
\alias{\%*\%,ndenseMatrix,matrix-method}
\alias{\%*\%,ndenseMatrix,ndenseMatrix-method}
\alias{\%*\%,dgeMatrix,dgeMatrix-method}
\alias{\%*\%,dgeMatrix,matrix-method}
\alias{\%*\%,dgeMatrix,numLike-method}
\alias{\%*\%,matrix,dgeMatrix-method}
\alias{\%*\%,numLike,dgeMatrix-method}
\alias{\%*\%,numLike,CsparseMatrix-method}
\alias{\%*\%,CsparseMatrix,CsparseMatrix-method}
\alias{\%*\%,CsparseMatrix,ddenseMatrix-method}
\alias{\%*\%,CsparseMatrix,matrix-method}
\alias{\%*\%,CsparseMatrix,numLike-method}
\alias{\%*\%,ddenseMatrix,CsparseMatrix-method}
\alias{\%*\%,matrix,CsparseMatrix-method}
\alias{\%*\%,matrix,sparseMatrix-method}
\alias{\%*\%,sparseMatrix,matrix-method}
\alias{\%*\%,ddenseMatrix,ddenseMatrix-method}
\alias{\%*\%,dgeMatrix,diagonalMatrix-method}
\alias{\%*\%,matrix,diagonalMatrix-method}
\alias{\%*\%,diagonalMatrix,dgeMatrix-method}
\alias{\%*\%,diagonalMatrix,diagonalMatrix-method}
\alias{\%*\%,diagonalMatrix,matrix-method}
\alias{\%*\%,diagonalMatrix,CsparseMatrix-method}
\alias{\%*\%,diagonalMatrix,sparseMatrix-method}
\alias{\%*\%,CsparseMatrix,diagonalMatrix-method}
\alias{\%*\%,sparseMatrix,diagonalMatrix-method}
\alias{\%*\%,denseMatrix,diagonalMatrix-method}
\alias{\%*\%,diagonalMatrix,denseMatrix-method}
\alias{\%*\%,diagonalMatrix,lgeMatrix-method}
\alias{\%*\%,lgeMatrix,diagonalMatrix-method}
\alias{\%*\%,nMatrix,lMatrix-method}
\alias{\%*\%,lMatrix,nMatrix-method}
\alias{\%*\%,lMatrix,lMatrix-method}
\alias{\%*\%,nMatrix,nMatrix-method}
\alias{\%*\%,dMatrix,lMatrix-method}
\alias{\%*\%,dMatrix,nMatrix-method}
\alias{\%*\%,lMatrix,dMatrix-method}
\alias{\%*\%,nMatrix,dMatrix-method}
\alias{\%*\%,dMatrix,integer-method}
\alias{\%*\%,integer,dMatrix-method}
\alias{\%*\%,Matrix,numLike-method}
\alias{\%*\%,numLike,Matrix-method}
\alias{\%*\%,Matrix,ANY-method}
\alias{\%*\%,ANY,Matrix-method}
\alias{\%*\%,Matrix,matrix-method}
\alias{\%*\%,matrix,Matrix-method}
\alias{\%*\%,ddenseMatrix,dsparseMatrix-method}
\alias{\%*\%,dgeMatrix,dsparseMatrix-method}
\alias{\%*\%,dsparseMatrix,ddenseMatrix-method}
\alias{\%*\%,dsparseMatrix,dgeMatrix-method}
\alias{\%*\%,ddenseMatrix,dsyMatrix-method}
\alias{\%*\%,matrix,dsyMatrix-method}
\alias{\%*\%,dsyMatrix,dsyMatrix-method}
\alias{\%*\%,dspMatrix,ddenseMatrix-method}
\alias{\%*\%,dspMatrix,matrix-method}
\alias{\%*\%,dsyMatrix,ddenseMatrix-method}
\alias{\%*\%,dsyMatrix,matrix-method}
\alias{\%*\%,dtpMatrix,ddenseMatrix-method}
\alias{\%*\%,dgeMatrix,dtpMatrix-method}
\alias{\%*\%,dtpMatrix,matrix-method}
\alias{\%*\%,matrix,dtpMatrix-method}
\alias{\%*\%,dtrMatrix,dtrMatrix-method}
\alias{\%*\%,ddenseMatrix,dtrMatrix-method}
\alias{\%*\%,dtrMatrix,ddenseMatrix-method}
\alias{\%*\%,dtrMatrix,matrix-method}
\alias{\%*\%,matrix,dtrMatrix-method}
\alias{\%*\%,matrix,indMatrix-method}
\alias{\%*\%,indMatrix,matrix-method}
\alias{\%*\%,indMatrix,indMatrix-method}
\alias{\%*\%,Matrix,indMatrix-method}
\alias{\%*\%,indMatrix,Matrix-method}
\alias{\%*\%,lgCMatrix,lgCMatrix-method}
\alias{\%*\%,lsparseMatrix,lsparseMatrix-method}
\alias{\%*\%,lsparseMatrix,ldenseMatrix-method}
\alias{\%*\%,ldenseMatrix,lsparseMatrix-method}
\alias{\%*\%,ngCMatrix,ngCMatrix-method}
\alias{\%*\%,nsparseMatrix,nsparseMatrix-method}
\alias{\%*\%,nsparseMatrix,ndenseMatrix-method}
\alias{\%*\%,ndenseMatrix,nsparseMatrix-method}
\alias{\%*\%,matrix,pMatrix-method}
\alias{\%*\%,pMatrix,matrix-method}
\alias{\%*\%,pMatrix,pMatrix-method}
\alias{\%*\%,Matrix,pMatrix-method}
\alias{\%*\%,pMatrix,Matrix-method}
\alias{\%*\%,mMatrix,sparseVector-method}
\alias{\%*\%,sparseVector,mMatrix-method}
\alias{\%*\%,sparseVector,sparseVector-method}
\alias{\%*\%,numLike,sparseVector-method}
\alias{\%*\%,sparseVector,numLike-method}
\alias{\%*\%,TsparseMatrix,ANY-method}
\alias{\%*\%,ANY,TsparseMatrix-method}
\alias{\%*\%,TsparseMatrix,Matrix-method}
\alias{\%*\%,Matrix,TsparseMatrix-method}
\alias{\%*\%,TsparseMatrix,TsparseMatrix-method}
%%-- crossprod -------------------------------
\alias{crossprod-methods}
\alias{crossprod,dgeMatrix,dgeMatrix-method}
\alias{crossprod,dgeMatrix,matrix-method}
\alias{crossprod,dgeMatrix,numLike-method}
\alias{crossprod,dgeMatrix,missing-method}
\alias{crossprod,matrix,dgeMatrix-method}
\alias{crossprod,numLike,dgeMatrix-method}
\alias{crossprod,numLike,CsparseMatrix-method}
\alias{crossprod,CsparseMatrix,missing-method}
\alias{crossprod,CsparseMatrix,CsparseMatrix-method}
\alias{crossprod,CsparseMatrix,ddenseMatrix-method}
\alias{crossprod,CsparseMatrix,matrix-method}
\alias{crossprod,CsparseMatrix,numLike-method}
\alias{crossprod,ddenseMatrix,CsparseMatrix-method}
\alias{crossprod,matrix,CsparseMatrix-method}
\alias{crossprod,ddenseMatrix,missing-method}
\alias{crossprod,ddenseMatrix,dgCMatrix-method}
\alias{crossprod,dgCMatrix,dgeMatrix-method}
\alias{crossprod,CsparseMatrix,diagonalMatrix-method}
\alias{crossprod,diagonalMatrix,CsparseMatrix-method}
\alias{crossprod,diagonalMatrix,dgeMatrix-method}
\alias{crossprod,diagonalMatrix,diagonalMatrix-method}
\alias{crossprod,diagonalMatrix,lgeMatrix-method}
\alias{crossprod,diagonalMatrix,matrix-method}
\alias{crossprod,diagonalMatrix,missing-method}
\alias{crossprod,diagonalMatrix,sparseMatrix-method}
\alias{crossprod,sparseMatrix,diagonalMatrix-method}
\alias{crossprod,ANY,Matrix-method}
\alias{crossprod,Matrix,numLike-method}
\alias{crossprod,numLike,Matrix-method}
\alias{crossprod,Matrix,ANY-method}
\alias{crossprod,Matrix,missing-method}
\alias{crossprod,Matrix,Matrix-method}
\alias{crossprod,Matrix,matrix-method}
\alias{crossprod,matrix,Matrix-method}
\alias{crossprod,matrix,diagonalMatrix-method}
\alias{crossprod,ddenseMatrix,dsparseMatrix-method}
\alias{crossprod,dgeMatrix,dsparseMatrix-method}
\alias{crossprod,dsparseMatrix,ddenseMatrix-method}
\alias{crossprod,dsparseMatrix,dgeMatrix-method}
\alias{crossprod,dtpMatrix,ddenseMatrix-method}
\alias{crossprod,dtpMatrix,matrix-method}
\alias{crossprod,dtrMatrix,ddenseMatrix-method}
\alias{crossprod,dtrMatrix,dtrMatrix-method}
\alias{crossprod,dtrMatrix,matrix-method}
\alias{crossprod,dtrMatrix,missing-method}
\alias{crossprod,indMatrix,matrix-method}
\alias{crossprod,indMatrix,Matrix-method}
\alias{crossprod,indMatrix,indMatrix-method}
\alias{crossprod,indMatrix,missing-method}
\alias{crossprod,lgCMatrix,missing-method}
\alias{crossprod,lgTMatrix,missing-method}
\alias{crossprod,lsparseMatrix-method}
\alias{crossprod,lsparseMatrix,missing-method}
\alias{crossprod,lsparseMatrix,lsparseMatrix-method}
\alias{crossprod,lsparseMatrix,ldenseMatrix-method}
\alias{crossprod,ldenseMatrix,lsparseMatrix-method}
\alias{crossprod,ngCMatrix,missing-method}
\alias{crossprod,ngTMatrix,missing-method}
\alias{crossprod,nsparseMatrix-method}
\alias{crossprod,nsparseMatrix,missing-method}
\alias{crossprod,nsparseMatrix,nsparseMatrix-method}
\alias{crossprod,nsparseMatrix,ndenseMatrix-method}
\alias{crossprod,ndenseMatrix,nsparseMatrix-method}
\alias{crossprod,pMatrix,matrix-method}
\alias{crossprod,pMatrix,Matrix-method}
\alias{crossprod,pMatrix,pMatrix-method}
\alias{crossprod,pMatrix,missing-method}
\alias{crossprod,mMatrix,sparseVector-method}
\alias{crossprod,sparseVector,mMatrix-method}
\alias{crossprod,sparseVector,missing-method}
\alias{crossprod,sparseVector,sparseVector-method}
\alias{crossprod,numLike,sparseVector-method}
\alias{crossprod,sparseVector,numLike-method}
\alias{crossprod,TsparseMatrix,ANY-method}
\alias{crossprod,ANY,TsparseMatrix-method}
\alias{crossprod,TsparseMatrix,Matrix-method}
\alias{crossprod,Matrix,TsparseMatrix-method}
\alias{crossprod,TsparseMatrix,TsparseMatrix-method}
\alias{crossprod,TsparseMatrix,missing-method}
\alias{crossprod,symmetricMatrix,Matrix-method}
\alias{crossprod,symmetricMatrix,missing-method}
\alias{crossprod,symmetricMatrix,ANY-method}
\alias{crossprod,ddenseMatrix,ddenseMatrix-method}
\alias{crossprod,ddenseMatrix,ldenseMatrix-method}
\alias{crossprod,ddenseMatrix,matrix-method}
\alias{crossprod,ddenseMatrix,ndenseMatrix-method}
\alias{crossprod,ldenseMatrix,ddenseMatrix-method}
\alias{crossprod,ldenseMatrix,ldenseMatrix-method}
\alias{crossprod,ldenseMatrix,matrix-method}
\alias{crossprod,ldenseMatrix,missing-method}
\alias{crossprod,ldenseMatrix,ndenseMatrix-method}
\alias{crossprod,ndenseMatrix,ddenseMatrix-method}
\alias{crossprod,ndenseMatrix,ldenseMatrix-method}
\alias{crossprod,ndenseMatrix,matrix-method}
\alias{crossprod,ndenseMatrix,missing-method}
\alias{crossprod,ndenseMatrix,ndenseMatrix-method}
%%-- tcrossprod ------------------------------
\alias{tcrossprod-methods}
\alias{tcrossprod,dgeMatrix,missing-method}
\alias{tcrossprod,dgeMatrix,dgeMatrix-method}
\alias{tcrossprod,dgeMatrix,matrix-method}
\alias{tcrossprod,dgeMatrix,numLike-method}
\alias{tcrossprod,matrix,dgeMatrix-method}
\alias{tcrossprod,numLike,dgeMatrix-method}
\alias{tcrossprod,CsparseMatrix,ddenseMatrix-method}
\alias{tcrossprod,CsparseMatrix,matrix-method}
\alias{tcrossprod,CsparseMatrix,numLike-method}
\alias{tcrossprod,ddenseMatrix,CsparseMatrix-method}
\alias{tcrossprod,matrix,CsparseMatrix-method}
\alias{tcrossprod,numLike,CsparseMatrix-method}
\alias{tcrossprod,CsparseMatrix,CsparseMatrix-method}
\alias{tcrossprod,CsparseMatrix,missing-method}
\alias{tcrossprod,ddenseMatrix,missing-method}
\alias{tcrossprod,CsparseMatrix,diagonalMatrix-method}
\alias{tcrossprod,dgeMatrix,diagonalMatrix-method}
\alias{tcrossprod,diagonalMatrix,CsparseMatrix-method}
\alias{tcrossprod,diagonalMatrix,diagonalMatrix-method}
\alias{tcrossprod,diagonalMatrix,matrix-method}
\alias{tcrossprod,diagonalMatrix,missing-method}
\alias{tcrossprod,diagonalMatrix,sparseMatrix-method}
\alias{tcrossprod,lgeMatrix,diagonalMatrix-method}
\alias{tcrossprod,matrix,diagonalMatrix-method}
\alias{tcrossprod,sparseMatrix,diagonalMatrix-method}
\alias{tcrossprod,ANY,Matrix-method}
\alias{tcrossprod,numLike,Matrix-method}
\alias{tcrossprod,Matrix,numLike-method}
\alias{tcrossprod,Matrix,ANY-method}
\alias{tcrossprod,Matrix,missing-method}
\alias{tcrossprod,Matrix,Matrix-method}
\alias{tcrossprod,Matrix,matrix-method}
\alias{tcrossprod,matrix,Matrix-method}
\alias{tcrossprod,ddenseMatrix,dtrMatrix-method}
\alias{tcrossprod,dtrMatrix,dtrMatrix-method}
\alias{tcrossprod,matrix,dtrMatrix-method}
\alias{tcrossprod,matrix,indMatrix-method}
\alias{tcrossprod,Matrix,indMatrix-method}
\alias{tcrossprod,indMatrix,indMatrix-method}
\alias{tcrossprod,indMatrix,missing-method}
\alias{tcrossprod,lgCMatrix,missing-method}
\alias{tcrossprod,lgTMatrix,missing-method}
\alias{tcrossprod,lsparseMatrix,missing-method}
\alias{tcrossprod,ngCMatrix,missing-method}
\alias{tcrossprod,ngTMatrix,missing-method}
\alias{tcrossprod,nsparseMatrix,missing-method}
\alias{tcrossprod,matrix,pMatrix-method}
\alias{tcrossprod,Matrix,pMatrix-method}
\alias{tcrossprod,pMatrix,pMatrix-method}
\alias{tcrossprod,pMatrix,missing-method}
\alias{tcrossprod,mMatrix,sparseVector-method}
\alias{tcrossprod,sparseMatrix,sparseVector-method}
\alias{tcrossprod,sparseVector,mMatrix-method}
\alias{tcrossprod,sparseVector,missing-method}
\alias{tcrossprod,sparseVector,sparseMatrix-method}
\alias{tcrossprod,sparseVector,sparseVector-method}
\alias{tcrossprod,numLike,sparseVector-method}
\alias{tcrossprod,sparseVector,numLike-method}
\alias{tcrossprod,TsparseMatrix,ANY-method}
\alias{tcrossprod,ANY,TsparseMatrix-method}
\alias{tcrossprod,TsparseMatrix,Matrix-method}
\alias{tcrossprod,Matrix,TsparseMatrix-method}
\alias{tcrossprod,TsparseMatrix,TsparseMatrix-method}
\alias{tcrossprod,TsparseMatrix,missing-method}
\alias{tcrossprod,Matrix,symmetricMatrix-method}
\alias{tcrossprod,ANY,symmetricMatrix-method}
\alias{tcrossprod,ddenseMatrix,ddenseMatrix-method}
\alias{tcrossprod,ddenseMatrix,ldenseMatrix-method}
\alias{tcrossprod,ddenseMatrix,matrix-method}
\alias{tcrossprod,ddenseMatrix,ndenseMatrix-method}
\alias{tcrossprod,ldenseMatrix,ddenseMatrix-method}
\alias{tcrossprod,ldenseMatrix,ldenseMatrix-method}
\alias{tcrossprod,ldenseMatrix,matrix-method}
\alias{tcrossprod,ldenseMatrix,missing-method}
\alias{tcrossprod,ldenseMatrix,ndenseMatrix-method}
\alias{tcrossprod,ndenseMatrix,ddenseMatrix-method}
\alias{tcrossprod,ndenseMatrix,ldenseMatrix-method}
\alias{tcrossprod,ndenseMatrix,matrix-method}
\alias{tcrossprod,ndenseMatrix,missing-method}
\alias{tcrossprod,ndenseMatrix,ndenseMatrix-method}
% (special class methods, added only for speedup:)
\alias{tcrossprod,ddenseMatrix,dsCMatrix-method}
\alias{tcrossprod,ddenseMatrix,lsCMatrix-method}
\alias{tcrossprod,ddenseMatrix,nsCMatrix-method}
\alias{tcrossprod,matrix,dsCMatrix-method}
\alias{tcrossprod,matrix,lsCMatrix-method}
\alias{tcrossprod,matrix,nsCMatrix-method}
%%---------------
%- these three are needed for R CMD check:
\alias{\%*\%}
\alias{crossprod}
\alias{tcrossprod}
%% no longer (don't want to mask 'base' functionality):
% \alias{tcrossprod,matrix,missing-method}
% \alias{tcrossprod,numeric,missing-method}
\description{
  The basic matrix product, \code{\%*\%} is implemented for all our
  \code{\linkS4class{Matrix}} and also for
  \code{\linkS4class{sparseVector}} classes, fully analogously to \R's
  base \code{matrix} and vector objects.

  The functions \code{\link{crossprod}} and \code{\link{tcrossprod}} are
  matrix products or \dQuote{cross products}, ideally implemented
  efficiently without computing \code{\link{t}(.)}'s unnecessarily.
  They also return \code{\linkS4class{symmetricMatrix}} classed
  matrices when easily detectable, e.g., in \code{crossprod(m)}, the one
  argument case.

  \code{tcrossprod()} takes the cross-product of the transpose of a matrix.
  \code{tcrossprod(x)} is formally equivalent to, but faster than, the
  call \code{x \%*\% t(x)}, and so is \code{tcrossprod(x, y)} instead of
  \code{x \%*\% t(y)}.

  \emph{Boolean} matrix products are computed via either
  \code{\link{\%&\%}} or \code{boolArith = TRUE}.
}
\usage{
\S4method{\%*\%}{CsparseMatrix,diagonalMatrix}(x, y)

\S4method{crossprod}{dgeMatrix,missing}(x, y = NULL, boolArith = NA, \dots)
\S4method{crossprod}{CsparseMatrix,diagonalMatrix}(x, y = NULL, boolArith = NA, \dots)
       ## .... and for many more signatures

\S4method{tcrossprod}{CsparseMatrix,ddenseMatrix}(x, y = NULL, boolArith = NA, \dots)
\S4method{tcrossprod}{TsparseMatrix,missing}(x, y = NULL, boolArith = NA, \dots)
       ## .... and for many more signatures
}
\arguments{
  \item{x}{a matrix-like object}
  \item{y}{a matrix-like object, or for \code{[t]crossprod()}
    \code{NULL} (by default); the latter case is formally equivalent to
    \code{y = x}.}
  \item{boolArith}{\code{\link{logical}}, i.e., \code{NA}, \code{TRUE},
    or \code{FALSE}.  If true the result is (coerced to) a pattern
    matrix, i.e., \code{"\linkS4class{nMatrix}"}, unless there are
    \code{NA} entries and the result will be a
    \code{"\linkS4class{lMatrix}"}.  If false the result is (coerced to)
    numeric.  When \code{NA}, currently the default, the
    result is a pattern matrix when \code{x} and \code{y} are
    \code{"\linkS4class{nsparseMatrix}"} and numeric otherwise.}
  \item{\dots}{potentially more arguments passed to and from methods.}
}
\details{
  For some classes in the \code{Matrix} package, such as
  \code{\linkS4class{dgCMatrix}}, it is much faster to calculate the
  cross-product of the transpose directly instead of calculating the
  transpose first and then its cross-product.

  \code{boolArith = TRUE} for regular (\dQuote{non cross}) matrix
  products, \code{\%*\%} cannot be specified.  Instead, we provide the
  \code{\link{\%&\%}} operator for \emph{boolean} matrix products.
}
\note{
  \code{boolArith = TRUE}, \code{FALSE} or \code{NA} has been newly
  introduced for \pkg{Matrix} 1.2.0 (March 2015).  Its implementation
  may be incomplete and partly missing.  Please report such omissions if
  detected!

  Currently, \code{boolArith = TRUE} is implemnted via
  \code{\linkS4class{CsparseMatrix}} coercions which may be quite
  inefficient for dense matrices.  Contributions for efficiency
  improvements are welcome.
}
\value{
  A \code{\linkS4class{Matrix}} object, in the one argument case
  of an appropriate symmetric matrix class.
}
\section{Methods}{
  \describe{
    \item{\%*\%}{\code{signature(x = "dgeMatrix", y = "dgeMatrix")}:
      Matrix multiplication; ditto for several other signature
      combinations, see \code{showMethods("\%*\%", class = "dgeMatrix")}.}
    \item{\%*\%}{\code{signature(x = "dtrMatrix", y = "matrix")} and other
      signatures (use \code{showMethods("\%*\%", class="dtrMatrix")}):
      matrix multiplication.  Multiplication of (matching) triangular
      matrices now should remain triangular (in the sense of class
      \linkS4class{triangularMatrix}).}
    \item{crossprod}{\code{signature(x = "dgeMatrix", y = "dgeMatrix")}:
      ditto for several other signatures, use
      \code{showMethods("crossprod", class = "dgeMatrix")}, matrix
      crossproduct, an efficient version of \code{t(x) \%*\% y}.}
    \item{crossprod}{\code{signature(x = "CsparseMatrix", y = "missing")}
      returns \code{t(x) \%*\% x} as an \code{dsCMatrix} object.}
    \item{crossprod}{\code{signature(x = "TsparseMatrix", y = "missing")}
      returns \code{t(x) \%*\% x} as an \code{dsCMatrix} object.}

    \item{crossprod,tcrossprod}{\code{signature(x = "dtrMatrix", y =
      "matrix")} and other signatures, see \code{"\%*\%"} above.}
  }%{describe}
}
\seealso{\code{\link[base]{tcrossprod}} in \R's base,
  \code{\link{crossprod}} and \code{\link{\%*\%}}.}
\examples{
 ## A random sparse "incidence" matrix :
 m <- matrix(0, 400, 500)
 set.seed(12)
 m[runif(314, 0, length(m))] <- 1
 mm <- as(m, "dgCMatrix")
 object.size(m) / object.size(mm) # smaller by a factor of > 200

 ## tcrossprod() is very fast:
 system.time(tCmm <- tcrossprod(mm))# 0   (PIII, 933 MHz)
 system.time(cm <- crossprod(t(m))) # 0.16
 system.time(cm. <- tcrossprod(m))  # 0.02

 stopifnot(cm == as(tCmm, "matrix"))

 ## show sparse sub matrix
 tCmm[1:16, 1:30]
}
\keyword{methods}
\keyword{algebra}
