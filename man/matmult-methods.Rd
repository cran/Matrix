\name{matmult-methods}
\title{Matrix (Cross) Products (of Transpose)}
%
\docType{methods}
\keyword{algebra}
\keyword{arith}
\keyword{array}
%
\alias{\%*\%}
\alias{\%*\%-methods}
\alias{crossprod}
\alias{crossprod-methods}
\alias{tcrossprod}
\alias{tcrossprod-methods}
\alias{matmult-methods}
% %*%
\alias{\%*\%,ANY,Matrix-method}
\alias{\%*\%,ANY,sparseVector-method}
\alias{\%*\%,CsparseMatrix,CsparseMatrix-method}
\alias{\%*\%,CsparseMatrix,RsparseMatrix-method}
\alias{\%*\%,CsparseMatrix,TsparseMatrix-method}
\alias{\%*\%,CsparseMatrix,denseMatrix-method}
\alias{\%*\%,CsparseMatrix,diagonalMatrix-method}
\alias{\%*\%,CsparseMatrix,matrix-method}
\alias{\%*\%,CsparseMatrix,vector-method}
\alias{\%*\%,Matrix,ANY-method}
\alias{\%*\%,Matrix,indMatrix-method}
\alias{\%*\%,Matrix,pMatrix-method}
\alias{\%*\%,Matrix,sparseVector-method}
\alias{\%*\%,RsparseMatrix,CsparseMatrix-method}
\alias{\%*\%,RsparseMatrix,RsparseMatrix-method}
\alias{\%*\%,RsparseMatrix,TsparseMatrix-method}
\alias{\%*\%,RsparseMatrix,denseMatrix-method}
\alias{\%*\%,RsparseMatrix,diagonalMatrix-method}
\alias{\%*\%,RsparseMatrix,matrix-method}
\alias{\%*\%,RsparseMatrix,vector-method}
\alias{\%*\%,TsparseMatrix,CsparseMatrix-method}
\alias{\%*\%,TsparseMatrix,RsparseMatrix-method}
\alias{\%*\%,TsparseMatrix,TsparseMatrix-method}
\alias{\%*\%,TsparseMatrix,denseMatrix-method}
\alias{\%*\%,TsparseMatrix,diagonalMatrix-method}
\alias{\%*\%,TsparseMatrix,matrix-method}
\alias{\%*\%,TsparseMatrix,vector-method}
\alias{\%*\%,denseMatrix,CsparseMatrix-method}
\alias{\%*\%,denseMatrix,RsparseMatrix-method}
\alias{\%*\%,denseMatrix,TsparseMatrix-method}
\alias{\%*\%,denseMatrix,denseMatrix-method}
\alias{\%*\%,denseMatrix,diagonalMatrix-method}
\alias{\%*\%,denseMatrix,matrix-method}
\alias{\%*\%,denseMatrix,vector-method}
\alias{\%*\%,diagonalMatrix,CsparseMatrix-method}
\alias{\%*\%,diagonalMatrix,RsparseMatrix-method}
\alias{\%*\%,diagonalMatrix,TsparseMatrix-method}
\alias{\%*\%,diagonalMatrix,denseMatrix-method}
\alias{\%*\%,diagonalMatrix,diagonalMatrix-method}
\alias{\%*\%,diagonalMatrix,matrix-method}
\alias{\%*\%,diagonalMatrix,vector-method}
\alias{\%*\%,indMatrix,Matrix-method}
\alias{\%*\%,indMatrix,indMatrix-method}
\alias{\%*\%,indMatrix,matrix-method}
\alias{\%*\%,indMatrix,pMatrix-method}
\alias{\%*\%,indMatrix,vector-method}
\alias{\%*\%,matrix,CsparseMatrix-method}
\alias{\%*\%,matrix,RsparseMatrix-method}
\alias{\%*\%,matrix,TsparseMatrix-method}
\alias{\%*\%,matrix,denseMatrix-method}
\alias{\%*\%,matrix,diagonalMatrix-method}
\alias{\%*\%,matrix,indMatrix-method}
\alias{\%*\%,matrix,pMatrix-method}
\alias{\%*\%,matrix,sparseVector-method}
\alias{\%*\%,pMatrix,Matrix-method}
\alias{\%*\%,pMatrix,indMatrix-method}
\alias{\%*\%,pMatrix,matrix-method}
\alias{\%*\%,pMatrix,pMatrix-method}
\alias{\%*\%,pMatrix,vector-method}
\alias{\%*\%,sparseVector,ANY-method}
\alias{\%*\%,sparseVector,Matrix-method}
\alias{\%*\%,sparseVector,matrix-method}
\alias{\%*\%,sparseVector,sparseVector-method}
\alias{\%*\%,sparseVector,vector-method}
\alias{\%*\%,vector,CsparseMatrix-method}
\alias{\%*\%,vector,RsparseMatrix-method}
\alias{\%*\%,vector,TsparseMatrix-method}
\alias{\%*\%,vector,denseMatrix-method}
\alias{\%*\%,vector,diagonalMatrix-method}
\alias{\%*\%,vector,indMatrix-method}
\alias{\%*\%,vector,pMatrix-method}
\alias{\%*\%,vector,sparseVector-method}
% crossprod
\alias{crossprod,ANY,Matrix-method}
\alias{crossprod,ANY,sparseVector-method}
\alias{crossprod,CsparseMatrix,CsparseMatrix-method}
\alias{crossprod,CsparseMatrix,RsparseMatrix-method}
\alias{crossprod,CsparseMatrix,TsparseMatrix-method}
\alias{crossprod,CsparseMatrix,denseMatrix-method}
\alias{crossprod,CsparseMatrix,diagonalMatrix-method}
\alias{crossprod,CsparseMatrix,matrix-method}
\alias{crossprod,CsparseMatrix,missing-method}
\alias{crossprod,CsparseMatrix,vector-method}
\alias{crossprod,Matrix,ANY-method}
\alias{crossprod,Matrix,indMatrix-method}
\alias{crossprod,Matrix,pMatrix-method}
\alias{crossprod,Matrix,sparseVector-method}
\alias{crossprod,RsparseMatrix,CsparseMatrix-method}
\alias{crossprod,RsparseMatrix,RsparseMatrix-method}
\alias{crossprod,RsparseMatrix,TsparseMatrix-method}
\alias{crossprod,RsparseMatrix,denseMatrix-method}
\alias{crossprod,RsparseMatrix,diagonalMatrix-method}
\alias{crossprod,RsparseMatrix,matrix-method}
\alias{crossprod,RsparseMatrix,missing-method}
\alias{crossprod,RsparseMatrix,vector-method}
\alias{crossprod,TsparseMatrix,CsparseMatrix-method}
\alias{crossprod,TsparseMatrix,RsparseMatrix-method}
\alias{crossprod,TsparseMatrix,TsparseMatrix-method}
\alias{crossprod,TsparseMatrix,denseMatrix-method}
\alias{crossprod,TsparseMatrix,diagonalMatrix-method}
\alias{crossprod,TsparseMatrix,matrix-method}
\alias{crossprod,TsparseMatrix,missing-method}
\alias{crossprod,TsparseMatrix,vector-method}
\alias{crossprod,denseMatrix,CsparseMatrix-method}
\alias{crossprod,denseMatrix,RsparseMatrix-method}
\alias{crossprod,denseMatrix,TsparseMatrix-method}
\alias{crossprod,denseMatrix,denseMatrix-method}
\alias{crossprod,denseMatrix,diagonalMatrix-method}
\alias{crossprod,denseMatrix,matrix-method}
\alias{crossprod,denseMatrix,missing-method}
\alias{crossprod,denseMatrix,vector-method}
\alias{crossprod,diagonalMatrix,CsparseMatrix-method}
\alias{crossprod,diagonalMatrix,RsparseMatrix-method}
\alias{crossprod,diagonalMatrix,TsparseMatrix-method}
\alias{crossprod,diagonalMatrix,denseMatrix-method}
\alias{crossprod,diagonalMatrix,diagonalMatrix-method}
\alias{crossprod,diagonalMatrix,matrix-method}
\alias{crossprod,diagonalMatrix,missing-method}
\alias{crossprod,diagonalMatrix,vector-method}
\alias{crossprod,indMatrix,Matrix-method}
\alias{crossprod,indMatrix,matrix-method}
\alias{crossprod,indMatrix,missing-method}
\alias{crossprod,indMatrix,vector-method}
\alias{crossprod,matrix,CsparseMatrix-method}
\alias{crossprod,matrix,RsparseMatrix-method}
\alias{crossprod,matrix,TsparseMatrix-method}
\alias{crossprod,matrix,denseMatrix-method}
\alias{crossprod,matrix,diagonalMatrix-method}
\alias{crossprod,matrix,indMatrix-method}
\alias{crossprod,matrix,pMatrix-method}
\alias{crossprod,matrix,sparseVector-method}
\alias{crossprod,pMatrix,missing-method}
\alias{crossprod,pMatrix,pMatrix-method}
\alias{crossprod,sparseVector,ANY-method}
\alias{crossprod,sparseVector,Matrix-method}
\alias{crossprod,sparseVector,matrix-method}
\alias{crossprod,sparseVector,missing-method}
\alias{crossprod,sparseVector,sparseVector-method}
\alias{crossprod,sparseVector,vector-method}
\alias{crossprod,vector,CsparseMatrix-method}
\alias{crossprod,vector,RsparseMatrix-method}
\alias{crossprod,vector,TsparseMatrix-method}
\alias{crossprod,vector,denseMatrix-method}
\alias{crossprod,vector,diagonalMatrix-method}
\alias{crossprod,vector,indMatrix-method}
\alias{crossprod,vector,pMatrix-method}
\alias{crossprod,vector,sparseVector-method}
% tcrossprod
\alias{tcrossprod,ANY,Matrix-method}
\alias{tcrossprod,ANY,sparseVector-method}
\alias{tcrossprod,CsparseMatrix,CsparseMatrix-method}
\alias{tcrossprod,CsparseMatrix,RsparseMatrix-method}
\alias{tcrossprod,CsparseMatrix,TsparseMatrix-method}
\alias{tcrossprod,CsparseMatrix,denseMatrix-method}
\alias{tcrossprod,CsparseMatrix,diagonalMatrix-method}
\alias{tcrossprod,CsparseMatrix,matrix-method}
\alias{tcrossprod,CsparseMatrix,missing-method}
\alias{tcrossprod,CsparseMatrix,vector-method}
\alias{tcrossprod,Matrix,ANY-method}
\alias{tcrossprod,Matrix,indMatrix-method}
\alias{tcrossprod,Matrix,sparseVector-method}
\alias{tcrossprod,RsparseMatrix,CsparseMatrix-method}
\alias{tcrossprod,RsparseMatrix,RsparseMatrix-method}
\alias{tcrossprod,RsparseMatrix,TsparseMatrix-method}
\alias{tcrossprod,RsparseMatrix,denseMatrix-method}
\alias{tcrossprod,RsparseMatrix,diagonalMatrix-method}
\alias{tcrossprod,RsparseMatrix,matrix-method}
\alias{tcrossprod,RsparseMatrix,missing-method}
\alias{tcrossprod,RsparseMatrix,vector-method}
\alias{tcrossprod,TsparseMatrix,CsparseMatrix-method}
\alias{tcrossprod,TsparseMatrix,RsparseMatrix-method}
\alias{tcrossprod,TsparseMatrix,TsparseMatrix-method}
\alias{tcrossprod,TsparseMatrix,denseMatrix-method}
\alias{tcrossprod,TsparseMatrix,diagonalMatrix-method}
\alias{tcrossprod,TsparseMatrix,matrix-method}
\alias{tcrossprod,TsparseMatrix,missing-method}
\alias{tcrossprod,TsparseMatrix,vector-method}
\alias{tcrossprod,denseMatrix,CsparseMatrix-method}
\alias{tcrossprod,denseMatrix,RsparseMatrix-method}
\alias{tcrossprod,denseMatrix,TsparseMatrix-method}
\alias{tcrossprod,denseMatrix,denseMatrix-method}
\alias{tcrossprod,denseMatrix,diagonalMatrix-method}
\alias{tcrossprod,denseMatrix,matrix-method}
\alias{tcrossprod,denseMatrix,missing-method}
\alias{tcrossprod,denseMatrix,vector-method}
\alias{tcrossprod,diagonalMatrix,CsparseMatrix-method}
\alias{tcrossprod,diagonalMatrix,RsparseMatrix-method}
\alias{tcrossprod,diagonalMatrix,TsparseMatrix-method}
\alias{tcrossprod,diagonalMatrix,denseMatrix-method}
\alias{tcrossprod,diagonalMatrix,diagonalMatrix-method}
\alias{tcrossprod,diagonalMatrix,matrix-method}
\alias{tcrossprod,diagonalMatrix,missing-method}
\alias{tcrossprod,diagonalMatrix,vector-method}
\alias{tcrossprod,indMatrix,Matrix-method}
\alias{tcrossprod,indMatrix,matrix-method}
\alias{tcrossprod,indMatrix,missing-method}
\alias{tcrossprod,indMatrix,vector-method}
\alias{tcrossprod,matrix,CsparseMatrix-method}
\alias{tcrossprod,matrix,RsparseMatrix-method}
\alias{tcrossprod,matrix,TsparseMatrix-method}
\alias{tcrossprod,matrix,denseMatrix-method}
\alias{tcrossprod,matrix,diagonalMatrix-method}
\alias{tcrossprod,matrix,indMatrix-method}
\alias{tcrossprod,matrix,sparseVector-method}
\alias{tcrossprod,pMatrix,Matrix-method}
\alias{tcrossprod,pMatrix,matrix-method}
\alias{tcrossprod,pMatrix,missing-method}
\alias{tcrossprod,pMatrix,pMatrix-method}
\alias{tcrossprod,pMatrix,vector-method}
\alias{tcrossprod,sparseVector,ANY-method}
\alias{tcrossprod,sparseVector,Matrix-method}
\alias{tcrossprod,sparseVector,matrix-method}
\alias{tcrossprod,sparseVector,missing-method}
\alias{tcrossprod,sparseVector,sparseVector-method}
\alias{tcrossprod,sparseVector,vector-method}
\alias{tcrossprod,vector,CsparseMatrix-method}
\alias{tcrossprod,vector,RsparseMatrix-method}
\alias{tcrossprod,vector,TsparseMatrix-method}
\alias{tcrossprod,vector,denseMatrix-method}
\alias{tcrossprod,vector,diagonalMatrix-method}
\alias{tcrossprod,vector,indMatrix-method}
\alias{tcrossprod,vector,sparseVector-method}
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

\S4method{crossprod}{CsparseMatrix,diagonalMatrix}(x, y = NULL, boolArith = NA, \dots)
       ## .... and for many more signatures

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
\dontshow{ % for R_DEFAULT_PACKAGES=NULL
library(stats, pos = "package:base", verbose = FALSE)
library(utils, pos = "package:base", verbose = FALSE)
}
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
