\name{qr-methods}
\docType{methods}
\alias{qr}
\alias{qrR}
\alias{qr-methods}
\alias{qr,dgCMatrix-method}
\alias{qr,sparseMatrix-method}
\title{QR Decomposition -- S4 Methods and Generic}
\description{
  The \code{"Matrix"} package provides methods for the QR decomposition
  of special classes of matrices.  There is a generic function which uses
  \code{\link[base]{qr}} as default, but methods defined in this package
  can take extra arguments.  In particular there is an option for
  determining a fill-reducing permutation of the columns of a sparse,
  rectangular matrix.
}
\usage{
qr(x, \dots)
qrR(qr, complete=FALSE, backPermute=TRUE, row.names = TRUE)
}
\arguments{
  \item{x}{a numeric or complex matrix whose QR decomposition is to be
    computed.  Logical matrices are coerced to numeric.}
  \item{qr}{a QR decomposition of the type computed by \code{qr}.}
  \item{complete}{logical indicating whether the \eqn{\bold{R}}
    matrix is to be completed by binding zero-value rows beneath the
    square upper triangle.}
  \item{backPermute}{logical indicating if the rows of the \eqn{\bold{R}}
    matrix should be back permuted such that \code{qrR()}'s result can
    be used directly to reconstruct the original matrix \eqn{\bold{X}}.}
  \item{row.names}{logical indicating if \code{\link{rownames}} should
    propagated to the result.}
  \item{\dots}{further arguments passed to or from other methods}
}
\section{Methods}{
  \describe{
    \item{x = "dgCMatrix"}{QR decomposition of a general sparse
      double-precision matrix with \code{nrow(x) >= ncol(x)}.  Returns
      an object of class \code{"\linkS4class{sparseQR}"}.}
    \item{x = "sparseMatrix"}{works via \code{"dgCMatrix"}.}
  }
}
\seealso{\code{\link[base]{qr}}; then, the class documentations,
  mainly \code{\linkS4class{sparseQR}}, and also
  \code{\linkS4class{dgCMatrix}}.
}
\examples{
%% FIXME: This is currently mixing example + regression tests
##------------- example of pivoting -- from base'  qraux.Rd -------------
X <- Matrix(cbind(int = 1,
                  b1=rep(1:0, each=3), b2=rep(0:1, each=3),
                  c1=rep(c(1,0,0), 2), c2=rep(c(0,1,0), 2), c3=rep(c(0,0,1),2)),
            sparse=TRUE)
rownames(X) <- paste0("r", seq_len(nrow(X)))
dnX <- dimnames(X)
X # is singular, columns "b2" and "c3" are "extra"
c(rankMatrix(X)) # = 4 (not 6)
##----- regular case ------------------------------------------
Xr <- X[ , -c(3,6)] # the "regular" (non-singular) version of X
stopifnot(rankMatrix(Xr) == ncol(Xr))
Y <- cbind(y <- setNames(1:6, paste0("y", 1:6)))
## regular case:
m <- as.matrix
qXr <- qr(  Xr)
qxr <- qr(m(Xr))
qcfXy <- qr.coef (qXr, y)
qcfXY <- qr.coef (qXr, Y)
stopifnot(
          all.equal(qr.coef(qxr, y), cf <- c(int=6, b1=-3, c1=-2, c2=-1), tol=1e-15)
,
          all.equal(qr.coef(qxr, Y), as.matrix(cf), tol=1e-15)
,
          all.equal(unname(qcfXy), unname(cf), tol=1e-15) || # FAIL names: ## FIXME_______
          all.equal(qcfXy, cf, tol=1e-15)
,
          all.equal(unname(m(qcfXY)), unname(m(cf)), tol=1e-15) || # FAIL dimnames: ## FIXME_______
          all.equal(m(qcfXY), m(cf), tol=1e-15)
,
          all.equal(y, qr.fitted(qxr, y), tol=2e-15)
,
          all.equal(y, qr.fitted(qXr, y), tol=2e-15)
,
          all.equal(m(qr.fitted(qXr, Y)), qr.fitted(qxr, Y), tol=1e-15)
,
          all.equal(  qr.resid (qXr, y),  qr.resid (qxr, y), tol=1e-15)
,
          all.equal(m(qr.resid (qXr, Y)), qr.resid (qxr, Y), tol=1e-15)
)

##----- singular case -----------------------------------------------
(qX <- qr(  X))
 qx <- qr(m(X))
# both @p and @q are non-trivial permutations
stopifnot(identical(dimnames(X), dnX))# some versions changed X's dimnames!

drop0(R. <- qr.R(qX), tol=1e-15) # columns *permuted*: c3 b1 ..
Q. <- qr.Q(qX)
qI <- sort.list(qX@q) # the inverse 'q' permutation
(X. <- drop0(Q. \%*\% R.[, qI], tol=1e-15))## just = X, incl. correct colnames
stopifnot(all(X - X.) < 8*.Machine$double.eps,
          ## qrR(.) returns R already "back permuted" (as with qI):
          identical(R.[, qI], qrR(qX)) )
##
## In this sense, classical qr.coef() is fine:
cfqx <- qr.coef(qx, y) # quite different from
nna <- !is.na(cfqx)
stopifnot(all.equal(unname(qr.fitted(qx,y)),
                    as.numeric(X[,nna] \%*\% cfqx[nna])))
## FIXME: do these make *any* sense? --- should give warnings !
qr.coef(qX, y)
qr.coef(qX, Y)
rm(m)
}
\keyword{methods}
\keyword{algebra}
\keyword{array}
