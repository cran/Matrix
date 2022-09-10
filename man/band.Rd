\name{band}
\title{Extract bands of a matrix}
\docType{methods}
%
\alias{band}
\alias{band-methods}
\alias{triu}
\alias{triu-methods}
\alias{tril}
\alias{tril-methods}
%
\alias{band,CsparseMatrix-method}
\alias{band,RsparseMatrix-method}
\alias{band,TsparseMatrix-method}
\alias{band,denseMatrix-method}
\alias{band,diagonalMatrix-method}
\alias{band,indMatrix-method}
\alias{band,matrix-method}
\alias{triu,CsparseMatrix-method}
\alias{triu,RsparseMatrix-method}
\alias{triu,TsparseMatrix-method}
\alias{triu,denseMatrix-method}
\alias{triu,diagonalMatrix-method}
\alias{triu,indMatrix-method}
\alias{triu,matrix-method}
\alias{tril,CsparseMatrix-method}
\alias{tril,RsparseMatrix-method}
\alias{tril,TsparseMatrix-method}
\alias{tril,denseMatrix-method}
\alias{tril,diagonalMatrix-method}
\alias{tril,indMatrix-method}
\alias{tril,matrix-method}
%
\description{
  Return the matrix obtained by setting to zero elements below a diagonal
  (\code{triu}), above a diagonal (\code{tril}), or outside of a general
  band (\code{band}).
}
\usage{
band(x, k1, k2, \dots)
triu(x, k = 0, \dots)
tril(x, k = 0, \dots)
}
\arguments{
  \item{x}{a matrix-like object}
  \item{k,k1,k2}{integers specifying the diagonals that are not set to
    zero.  These are interpreted relative to the main diagonal, which
    is \code{k=0}.  Positive and negative values of \code{k} indicate
    diagonals above and below the main diagonal, respectively.}
  \item{\dots}{optional arguments passed methods (currently unused
    by package \pkg{Matrix})}
}
\details{
  \code{triu(x, k)} is equivalent to \code{band(x, k, dim(x)[2])}.
  Similarly,
  \code{tril(x, k)} is equivalent to \code{band(x, -dim(x)[1], k)}.
}
\value{
  An object of a suitable matrix class, inheriting from
  \code{\linkS4class{triangularMatrix}} where appropriate.
  It inherits from \code{\linkS4class{sparseMatrix}} if
  and only if \code{x} does.
}
\section{Methods}{
  \describe{
    \item{x = "CsparseMatrix"}{method for compressed, sparse,
      column-oriented matrices.}
    \item{x = "RsparseMatrix"}{method for compressed, sparse,
      row-oriented matrices.}
    \item{x = "TsparseMatrix"}{method for sparse matrices in
      triplet format.}
    \item{x = "diagonalMatrix"}{method for diagonal matrices.}
    \item{x = "denseMatrix"}{method for dense matrices in
      packed or unpacked format.}
    \item{x = "matrix"}{method for traditional matrices
      of implicit class \code{\link{matrix}}.}
  }
}
\seealso{
  \code{\link{bandSparse}} for the \emph{construction} of a
  banded sparse matrix directly from its non-zero diagonals.
}
\examples{
## A random sparse matrix :
set.seed(7)
m <- matrix(0, 5, 5)
m[sample(length(m), size = 14)] <- rep(1:9, length=14)
(mm <- as(m, "CsparseMatrix"))

tril(mm)        # lower triangle
tril(mm, -1)    # strict lower triangle
triu(mm,  1)    # strict upper triangle
band(mm, -1, 2) # general band
(m5 <- Matrix(rnorm(25), ncol = 5))
tril(m5)        # lower triangle
tril(m5, -1)    # strict lower triangle
triu(m5, 1)     # strict upper triangle
band(m5, -1, 2) # general band
(m65 <- Matrix(rnorm(30), ncol = 5))  # not square
triu(m65)       # result not "dtrMatrix" unless square
(sm5 <- crossprod(m65)) # symmetric
   band(sm5, -1, 1)# "dsyMatrix": symmetric band preserves symmetry property
as(band(sm5, -1, 1), "sparseMatrix")# often preferable
(sm <- round(crossprod(triu(mm/2)))) # sparse symmetric ("dsC*")
band(sm, -1,1) # remains "dsC", *however*
band(sm, -2,1) # -> "dgC"
%% Problem is the underlying  cholmod_band()  which really does symmetric
%% banding *only* *if* the matrix is cholmod-symmetric i.e. 'stype != 0'
%%
\dontshow{ ## this uses special methods
(x.x <- crossprod(mm))
tril(x.x)
xx <- tril(x.x) + triu(x.x, 1) ## the same as x.x (but stored differently):
txx <- t(as(xx, "symmetricMatrix"))
stopifnot(identical(triu(x.x), t(tril(x.x))),
	  identical(class(x.x), class(txx)),
	  identical(as(x.x, "generalMatrix"), as(txx, "generalMatrix")))
}
}
\keyword{methods}
\keyword{algebra}
