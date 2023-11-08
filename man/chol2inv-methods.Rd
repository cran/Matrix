\name{chol2inv-methods}
\title{Inverse from Cholesky Factor}
%
\docType{methods}
\keyword{algebra}
\keyword{array}
\keyword{methods}
%
\alias{chol2inv}
\alias{chol2inv-methods}
%
\alias{chol2inv,ANY-method}
\alias{chol2inv,ddiMatrix-method}
\alias{chol2inv,diagonalMatrix-method}
\alias{chol2inv,dtCMatrix-method}
\alias{chol2inv,dtRMatrix-method}
\alias{chol2inv,dtTMatrix-method}
\alias{chol2inv,dtrMatrix-method}
\alias{chol2inv,dtpMatrix-method}
\alias{chol2inv,generalMatrix-method}
\alias{chol2inv,symmetricMatrix-method}
\alias{chol2inv,triangularMatrix-method}
%
\description{
  Given \emph{formally} upper and lower triangular matrices
  \eqn{U} and \eqn{L}, compute \eqn{(U' U)^{-1}}{(U' U)^(-1)}
  and \eqn{(L L')^{-1}}{(L L')^(-1)}, respectively.

  This function can be seen as way to compute the inverse of a
  symmetric positive definite matrix given its Cholesky factor.
  Equivalently, it can be seen as a way to compute
  \eqn{(X' X)^{-1}}{(X' X)^(-1)} given the \eqn{R} part of the
  QR factorization of \eqn{X}, if \eqn{R} is constrained to have
  positive diagonal entries.
}
\usage{
chol2inv(x, \dots)
\S4method{chol2inv}{dtrMatrix}(x, \dots)
\S4method{chol2inv}{dtCMatrix}(x, \dots)
\S4method{chol2inv}{generalMatrix}(x, uplo = "U", \dots)
}
\arguments{
  \item{x}{a square matrix or \code{\linkS4class{Matrix}},
    typically the result of a call to \code{\link{chol}}.
    If \code{x} is square but not (formally) triangular,
    then only the upper or lower triangle is considered,
    depending on optional argument \code{uplo} if \code{x}
    is a \code{\linkS4class{Matrix}}.}
  \item{uplo}{a string, either \code{"U"} or \code{"L"},
    indicating which triangle of \code{x} contains the
    Cholesky factor.  The default is \code{"U"}, to be
    consistent with \code{\link[base]{chol2inv}} from \pkg{base}.}
  \item{\dots}{further arguments passed to or from methods.}
}
\value{
  A matrix, \code{\linkS4class{symmetricMatrix}},
  or \code{\linkS4class{diagonalMatrix}} representing
  the inverse of the positive definite matrix whose
  Cholesky factor is \code{x}.
  The result is a traditional matrix if \code{x} is a
  traditional matrix, dense if \code{x} is dense, and
  sparse if \code{x} is sparse.
}
\seealso{
  The default method from \pkg{base}, \code{\link[base]{chol2inv}},
  called for traditional matrices \code{x}.

  Generic function \code{\link{chol}}, for computing the upper
  triangular Cholesky factor \eqn{L'} of a symmetric positive
  semidefinite matrix.

  Generic function \code{\link{solve}}, for solving linear systems
  and (as a corollary) for computing inverses more generally.
}
\examples{
(A <- Matrix(cbind(c(1, 1, 1), c(1, 2, 4), c(1, 4, 16))))
(R <- chol(A))
(L <- t(R))
(R2i <- chol2inv(R))
(L2i <- chol2inv(R))
stopifnot(exprs = {
    all.equal(R2i, tcrossprod(solve(R)))
    all.equal(L2i,  crossprod(solve(L)))
    all.equal(as(R2i \%*\% A, "matrix"), diag(3L)) # the identity 
    all.equal(as(L2i \%*\% A, "matrix"), diag(3L)) # ditto
})
}
