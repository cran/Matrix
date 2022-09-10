\name{matrix-products}
\title{Matrix (Cross) Products (of Transpose)}
\docType{methods}
%
\alias{\%*\%-methods}
\alias{crossprod-methods}
\alias{tcrossprod-methods}
% and for R CMD check:
\alias{\%*\%}
\alias{crossprod}
\alias{tcrossprod}
% `%*%`
\alias{\%*\%,ANY,Matrix-method}
\alias{\%*\%,ANY,TsparseMatrix-method}
\alias{\%*\%,CsparseMatrix,CsparseMatrix-method}
\alias{\%*\%,CsparseMatrix,ddenseMatrix-method}
\alias{\%*\%,CsparseMatrix,diagonalMatrix-method}
\alias{\%*\%,CsparseMatrix,matrix-method}
\alias{\%*\%,CsparseMatrix,numLike-method}
\alias{\%*\%,Matrix,ANY-method}
\alias{\%*\%,Matrix,TsparseMatrix-method}
\alias{\%*\%,Matrix,indMatrix-method}
\alias{\%*\%,Matrix,matrix-method}
\alias{\%*\%,Matrix,numLike-method}
\alias{\%*\%,Matrix,pMatrix-method}
\alias{\%*\%,RsparseMatrix,mMatrix-method}
\alias{\%*\%,TsparseMatrix,ANY-method}
\alias{\%*\%,TsparseMatrix,Matrix-method}
\alias{\%*\%,TsparseMatrix,TsparseMatrix-method}
\alias{\%*\%,dMatrix,lMatrix-method}
\alias{\%*\%,dMatrix,nMatrix-method}
\alias{\%*\%,ddenseMatrix,CsparseMatrix-method}
\alias{\%*\%,ddenseMatrix,ddenseMatrix-method}
\alias{\%*\%,ddenseMatrix,dsyMatrix-method}
\alias{\%*\%,ddenseMatrix,dtrMatrix-method}
\alias{\%*\%,ddenseMatrix,ldenseMatrix-method}
\alias{\%*\%,ddenseMatrix,matrix-method}
\alias{\%*\%,ddenseMatrix,ndenseMatrix-method}
\alias{\%*\%,denseMatrix,diagonalMatrix-method}
\alias{\%*\%,dgeMatrix,dgeMatrix-method}
\alias{\%*\%,dgeMatrix,diagonalMatrix-method}
\alias{\%*\%,dgeMatrix,dtpMatrix-method}
\alias{\%*\%,dgeMatrix,matrix-method}
\alias{\%*\%,diagonalMatrix,CsparseMatrix-method}
\alias{\%*\%,diagonalMatrix,denseMatrix-method}
\alias{\%*\%,diagonalMatrix,dgeMatrix-method}
\alias{\%*\%,diagonalMatrix,diagonalMatrix-method}
\alias{\%*\%,diagonalMatrix,lgeMatrix-method}
\alias{\%*\%,diagonalMatrix,matrix-method}
\alias{\%*\%,diagonalMatrix,sparseMatrix-method}
\alias{\%*\%,dspMatrix,ddenseMatrix-method}
\alias{\%*\%,dspMatrix,matrix-method}
\alias{\%*\%,dsyMatrix,ddenseMatrix-method}
\alias{\%*\%,dsyMatrix,dsyMatrix-method}
\alias{\%*\%,dsyMatrix,matrix-method}
\alias{\%*\%,dtpMatrix,ddenseMatrix-method}
\alias{\%*\%,dtpMatrix,matrix-method}
\alias{\%*\%,dtrMatrix,ddenseMatrix-method}
\alias{\%*\%,dtrMatrix,dtrMatrix-method}
\alias{\%*\%,dtrMatrix,matrix-method}
\alias{\%*\%,indMatrix,Matrix-method}
\alias{\%*\%,indMatrix,matrix-method}
\alias{\%*\%,lMatrix,dMatrix-method}
\alias{\%*\%,lMatrix,lMatrix-method}
\alias{\%*\%,lMatrix,nMatrix-method}
\alias{\%*\%,ldenseMatrix,ddenseMatrix-method}
\alias{\%*\%,ldenseMatrix,ldenseMatrix-method}
\alias{\%*\%,ldenseMatrix,lsparseMatrix-method}
\alias{\%*\%,ldenseMatrix,matrix-method}
\alias{\%*\%,ldenseMatrix,ndenseMatrix-method}
\alias{\%*\%,lgeMatrix,diagonalMatrix-method}
\alias{\%*\%,lsparseMatrix,ldenseMatrix-method}
\alias{\%*\%,lsparseMatrix,lsparseMatrix-method}
\alias{\%*\%,mMatrix,RsparseMatrix-method}
\alias{\%*\%,mMatrix,sparseVector-method}
\alias{\%*\%,matrix,CsparseMatrix-method}
\alias{\%*\%,matrix,Matrix-method}
\alias{\%*\%,matrix,ddenseMatrix-method}
\alias{\%*\%,matrix,dgeMatrix-method}
\alias{\%*\%,matrix,diagonalMatrix-method}
\alias{\%*\%,matrix,dsyMatrix-method}
\alias{\%*\%,matrix,dtpMatrix-method}
\alias{\%*\%,matrix,dtrMatrix-method}
\alias{\%*\%,matrix,indMatrix-method}
\alias{\%*\%,matrix,ldenseMatrix-method}
\alias{\%*\%,matrix,ndenseMatrix-method}
\alias{\%*\%,matrix,pMatrix-method}
\alias{\%*\%,matrix,sparseMatrix-method}
\alias{\%*\%,nMatrix,dMatrix-method}
\alias{\%*\%,nMatrix,lMatrix-method}
\alias{\%*\%,nMatrix,nMatrix-method}
\alias{\%*\%,ndenseMatrix,ddenseMatrix-method}
\alias{\%*\%,ndenseMatrix,ldenseMatrix-method}
\alias{\%*\%,ndenseMatrix,matrix-method}
\alias{\%*\%,ndenseMatrix,ndenseMatrix-method}
\alias{\%*\%,ndenseMatrix,nsparseMatrix-method}
\alias{\%*\%,nsparseMatrix,ndenseMatrix-method}
\alias{\%*\%,nsparseMatrix,nsparseMatrix-method}
\alias{\%*\%,numLike,CsparseMatrix-method}
\alias{\%*\%,numLike,Matrix-method}
\alias{\%*\%,numLike,sparseVector-method}
\alias{\%*\%,pMatrix,pMatrix-method}
\alias{\%*\%,sparseMatrix,diagonalMatrix-method}
\alias{\%*\%,sparseMatrix,matrix-method}
\alias{\%*\%,sparseVector,mMatrix-method}
\alias{\%*\%,sparseVector,numLike-method}
\alias{\%*\%,sparseVector,sparseVector-method}
% crossprod
\alias{crossprod,ANY,ANY-method}
\alias{crossprod,ANY,Matrix-method}
\alias{crossprod,ANY,RsparseMatrix-method}
\alias{crossprod,ANY,TsparseMatrix-method}
\alias{crossprod,CsparseMatrix,CsparseMatrix-method}
\alias{crossprod,CsparseMatrix,ddenseMatrix-method}
\alias{crossprod,CsparseMatrix,diagonalMatrix-method}
\alias{crossprod,CsparseMatrix,matrix-method}
\alias{crossprod,CsparseMatrix,missing-method}
\alias{crossprod,CsparseMatrix,numLike-method}
\alias{crossprod,Matrix,ANY-method}
\alias{crossprod,Matrix,Matrix-method}
\alias{crossprod,Matrix,TsparseMatrix-method}
\alias{crossprod,Matrix,matrix-method}
\alias{crossprod,Matrix,missing-method}
\alias{crossprod,Matrix,numLike-method}
\alias{crossprod,RsparseMatrix,ANY-method}
\alias{crossprod,RsparseMatrix,mMatrix-method}
\alias{crossprod,TsparseMatrix,ANY-method}
\alias{crossprod,TsparseMatrix,Matrix-method}
\alias{crossprod,TsparseMatrix,TsparseMatrix-method}
\alias{crossprod,TsparseMatrix,missing-method}
\alias{crossprod,ddenseMatrix,CsparseMatrix-method}
\alias{crossprod,ddenseMatrix,ddenseMatrix-method}
\alias{crossprod,ddenseMatrix,dgCMatrix-method}
\alias{crossprod,ddenseMatrix,dsparseMatrix-method}
\alias{crossprod,ddenseMatrix,ldenseMatrix-method}
\alias{crossprod,ddenseMatrix,matrix-method}
\alias{crossprod,ddenseMatrix,missing-method}
\alias{crossprod,ddenseMatrix,ndenseMatrix-method}
\alias{crossprod,dgCMatrix,dgeMatrix-method}
\alias{crossprod,dgeMatrix,dgeMatrix-method}
\alias{crossprod,dgeMatrix,matrix-method}
\alias{crossprod,dgeMatrix,missing-method}
\alias{crossprod,dgeMatrix,numLike-method}
\alias{crossprod,diagonalMatrix,CsparseMatrix-method}
\alias{crossprod,diagonalMatrix,dgeMatrix-method}
\alias{crossprod,diagonalMatrix,diagonalMatrix-method}
\alias{crossprod,diagonalMatrix,lgeMatrix-method}
\alias{crossprod,diagonalMatrix,matrix-method}
\alias{crossprod,diagonalMatrix,missing-method}
\alias{crossprod,diagonalMatrix,sparseMatrix-method}
\alias{crossprod,dsparseMatrix,ddenseMatrix-method}
\alias{crossprod,dsparseMatrix,dgeMatrix-method}
\alias{crossprod,dtpMatrix,ddenseMatrix-method}
\alias{crossprod,dtpMatrix,matrix-method}
\alias{crossprod,dtrMatrix,ddenseMatrix-method}
\alias{crossprod,dtrMatrix,dtrMatrix-method}
\alias{crossprod,dtrMatrix,matrix-method}
\alias{crossprod,indMatrix,Matrix-method}
\alias{crossprod,indMatrix,indMatrix-method}
\alias{crossprod,indMatrix,matrix-method}
\alias{crossprod,indMatrix,missing-method}
\alias{crossprod,ldenseMatrix,ddenseMatrix-method}
\alias{crossprod,ldenseMatrix,ldenseMatrix-method}
\alias{crossprod,ldenseMatrix,lsparseMatrix-method}
\alias{crossprod,ldenseMatrix,matrix-method}
\alias{crossprod,ldenseMatrix,missing-method}
\alias{crossprod,ldenseMatrix,ndenseMatrix-method}
\alias{crossprod,lsparseMatrix,ldenseMatrix-method}
\alias{crossprod,lsparseMatrix,lsparseMatrix-method}
\alias{crossprod,mMatrix,RsparseMatrix-method}
\alias{crossprod,mMatrix,sparseVector-method}
\alias{crossprod,matrix,CsparseMatrix-method}
\alias{crossprod,matrix,Matrix-method}
\alias{crossprod,matrix,dgeMatrix-method}
\alias{crossprod,matrix,diagonalMatrix-method}
\alias{crossprod,matrix,dtrMatrix-method}
\alias{crossprod,ndenseMatrix,ddenseMatrix-method}
\alias{crossprod,ndenseMatrix,ldenseMatrix-method}
\alias{crossprod,ndenseMatrix,matrix-method}
\alias{crossprod,ndenseMatrix,missing-method}
\alias{crossprod,ndenseMatrix,ndenseMatrix-method}
\alias{crossprod,ndenseMatrix,nsparseMatrix-method}
\alias{crossprod,nsparseMatrix,ndenseMatrix-method}
\alias{crossprod,nsparseMatrix,nsparseMatrix-method}
\alias{crossprod,numLike,CsparseMatrix-method}
\alias{crossprod,numLike,Matrix-method}
\alias{crossprod,numLike,dgeMatrix-method}
\alias{crossprod,numLike,sparseVector-method}
\alias{crossprod,pMatrix,Matrix-method}
\alias{crossprod,pMatrix,matrix-method}
\alias{crossprod,pMatrix,missing-method}
\alias{crossprod,pMatrix,pMatrix-method}
\alias{crossprod,sparseMatrix,diagonalMatrix-method}
\alias{crossprod,sparseVector,mMatrix-method}
\alias{crossprod,sparseVector,missing-method}
\alias{crossprod,sparseVector,numLike-method}
\alias{crossprod,sparseVector,sparseVector-method}
\alias{crossprod,symmetricMatrix,ANY-method}
\alias{crossprod,symmetricMatrix,Matrix-method}
\alias{crossprod,symmetricMatrix,missing-method}
% tcrossprod
\alias{tcrossprod,ANY,ANY-method}
\alias{tcrossprod,ANY,Matrix-method}
\alias{tcrossprod,ANY,RsparseMatrix-method}
\alias{tcrossprod,ANY,TsparseMatrix-method}
\alias{tcrossprod,ANY,symmetricMatrix-method}
\alias{tcrossprod,CsparseMatrix,CsparseMatrix-method}
\alias{tcrossprod,CsparseMatrix,ddenseMatrix-method}
\alias{tcrossprod,CsparseMatrix,diagonalMatrix-method}
\alias{tcrossprod,CsparseMatrix,matrix-method}
\alias{tcrossprod,CsparseMatrix,missing-method}
\alias{tcrossprod,CsparseMatrix,numLike-method}
\alias{tcrossprod,Matrix,ANY-method}
\alias{tcrossprod,Matrix,Matrix-method}
\alias{tcrossprod,Matrix,TsparseMatrix-method}
\alias{tcrossprod,Matrix,indMatrix-method}
\alias{tcrossprod,Matrix,matrix-method}
\alias{tcrossprod,Matrix,missing-method}
\alias{tcrossprod,Matrix,numLike-method}
\alias{tcrossprod,Matrix,symmetricMatrix-method}
\alias{tcrossprod,RsparseMatrix,ANY-method}
\alias{tcrossprod,RsparseMatrix,mMatrix-method}
\alias{tcrossprod,TsparseMatrix,ANY-method}
\alias{tcrossprod,TsparseMatrix,Matrix-method}
\alias{tcrossprod,TsparseMatrix,TsparseMatrix-method}
\alias{tcrossprod,TsparseMatrix,missing-method}
\alias{tcrossprod,ddenseMatrix,CsparseMatrix-method}
\alias{tcrossprod,ddenseMatrix,ddenseMatrix-method}
\alias{tcrossprod,ddenseMatrix,dsCMatrix-method}
\alias{tcrossprod,ddenseMatrix,dtrMatrix-method}
\alias{tcrossprod,ddenseMatrix,ldenseMatrix-method}
\alias{tcrossprod,ddenseMatrix,lsCMatrix-method}
\alias{tcrossprod,ddenseMatrix,matrix-method}
\alias{tcrossprod,ddenseMatrix,missing-method}
\alias{tcrossprod,ddenseMatrix,ndenseMatrix-method}
\alias{tcrossprod,ddenseMatrix,nsCMatrix-method}
\alias{tcrossprod,dgeMatrix,dgeMatrix-method}
\alias{tcrossprod,dgeMatrix,diagonalMatrix-method}
\alias{tcrossprod,dgeMatrix,matrix-method}
\alias{tcrossprod,dgeMatrix,missing-method}
\alias{tcrossprod,dgeMatrix,numLike-method}
\alias{tcrossprod,diagonalMatrix,CsparseMatrix-method}
\alias{tcrossprod,diagonalMatrix,diagonalMatrix-method}
\alias{tcrossprod,diagonalMatrix,missing-method}
\alias{tcrossprod,diagonalMatrix,sparseMatrix-method}
\alias{tcrossprod,dtrMatrix,dtrMatrix-method}
\alias{tcrossprod,indMatrix,indMatrix-method}
\alias{tcrossprod,indMatrix,missing-method}
\alias{tcrossprod,ldenseMatrix,ddenseMatrix-method}
\alias{tcrossprod,ldenseMatrix,ldenseMatrix-method}
\alias{tcrossprod,ldenseMatrix,matrix-method}
\alias{tcrossprod,ldenseMatrix,missing-method}
\alias{tcrossprod,ldenseMatrix,ndenseMatrix-method}
\alias{tcrossprod,lgeMatrix,diagonalMatrix-method}
\alias{tcrossprod,mMatrix,RsparseMatrix-method}
\alias{tcrossprod,mMatrix,sparseVector-method}
\alias{tcrossprod,matrix,CsparseMatrix-method}
\alias{tcrossprod,matrix,Matrix-method}
\alias{tcrossprod,matrix,dgeMatrix-method}
\alias{tcrossprod,matrix,diagonalMatrix-method}
\alias{tcrossprod,matrix,dsCMatrix-method}
\alias{tcrossprod,matrix,dtrMatrix-method}
\alias{tcrossprod,matrix,indMatrix-method}
\alias{tcrossprod,matrix,lsCMatrix-method}
\alias{tcrossprod,matrix,nsCMatrix-method}
\alias{tcrossprod,ndenseMatrix,ddenseMatrix-method}
\alias{tcrossprod,ndenseMatrix,ldenseMatrix-method}
\alias{tcrossprod,ndenseMatrix,matrix-method}
\alias{tcrossprod,ndenseMatrix,missing-method}
\alias{tcrossprod,ndenseMatrix,ndenseMatrix-method}
\alias{tcrossprod,numLike,CsparseMatrix-method}
\alias{tcrossprod,numLike,Matrix-method}
\alias{tcrossprod,numLike,dgeMatrix-method}
\alias{tcrossprod,numLike,sparseVector-method}
\alias{tcrossprod,pMatrix,missing-method}
\alias{tcrossprod,pMatrix,pMatrix-method}
\alias{tcrossprod,sparseMatrix,diagonalMatrix-method}
\alias{tcrossprod,sparseMatrix,sparseVector-method}
\alias{tcrossprod,sparseVector,mMatrix-method}
\alias{tcrossprod,sparseVector,missing-method}
\alias{tcrossprod,sparseVector,numLike-method}
\alias{tcrossprod,sparseVector,sparseMatrix-method}
\alias{tcrossprod,sparseVector,sparseVector-method}
% no longer, as not to mask base::tcrossprod
%\alias{tcrossprod,matrix,missing-method}
%\alias{tcrossprod,numeric,missing-method}
%
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
  has still not been tested extensively.  Notably the behaviour for
  sparse matrices with \code{x} slots containing extra zeros had not been
  documented previously, see the \code{\link{\%&\%}} help page.

  Currently, \code{boolArith = TRUE} is implemented via
  \code{\linkS4class{CsparseMatrix}} coercions which may be quite
  inefficient for dense matrices.  Contributions for efficiency
  improvements are welcome.
}
\value{
  A \code{\linkS4class{Matrix}} object, in the one argument case
  of an appropriate \emph{symmetric} matrix class, i.e., inheriting from
  \code{\linkS4class{symmetricMatrix}}.
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
\seealso{
  \code{\link[base]{tcrossprod}} in \R's base, and
  \code{\link{crossprod}} and \code{\link{\%*\%}}.
  \pkg{Matrix} package \code{\link{\%&\%}} for boolean matrix product
  methods.
}
\examples{
 ## A random sparse "incidence" matrix :
 m <- matrix(0, 400, 500)
 set.seed(12)
 m[runif(314, 0, length(m))] <- 1
 mm <- as(m, "CsparseMatrix")
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
