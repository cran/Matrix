\name{band}
\docType{methods}
\alias{band-methods}
\alias{tril-methods}
\alias{triu-methods}
\alias{band,CsparseMatrix-method}
\alias{tril,CsparseMatrix-method}
\alias{triu,CsparseMatrix-method}
\alias{band,RsparseMatrix-method}
\alias{tril,RsparseMatrix-method}
\alias{triu,RsparseMatrix-method}
\alias{band,TsparseMatrix-method}
\alias{tril,TsparseMatrix-method}
\alias{triu,TsparseMatrix-method}
\alias{tril,dsCMatrix-method}
\alias{tril,lsCMatrix-method}
\alias{triu,dsCMatrix-method}
\alias{triu,lsCMatrix-method}
%
\alias{band}
\alias{tril}
\alias{triu}
\title{Extract bands of a matrix}
\description{
  Returns a new matrix formed by extracting the lower triangle
  (\code{tril}) or the upper triangle (\code{triu}) or a general band
  relative to the diagonal (\code{band}), and setting other elements
  to zero.  The general forms of these functions include integer
  arguments to specify how many diagonal bands above or below the main
  diagonal are not set to zero.
}
\usage{
band(x, k1, k2, \dots)
tril(x, k = 0, \dots)
triu(x, k = 0, \dots)
}
\arguments{
  \item{x}{a matrix-like object}
  \item{k,k1,k2}{integers specifying the diagonal bands that will not
  be set to zero.  These are given relative to the main diagonal,
  which is \code{k=0}.  A negative value of \code{k} indicates a
  diagonal below the main diagonal and a positive value indicates a
  diagonal above the main diagonal.}
 \item{\dots}{Optional arguments used by specific methods. (None used at present.)}
}
\value{
  An object of an appropriate matrix class.  The class of the value of
  \code{tril} or \code{triu} inherits from
  \code{\linkS4class{triangularMatrix}} when appropriate.
}
\section{Methods}{
  \describe{
    \item{x = "CsparseMatrix"}{method for compressed, sparse,
      column-oriented matrices.}
    \item{x = "TsparseMatrix"}{method for sparse matrices in triplet format.}
    \item{x = "RsparseMatrix"}{method for compressed, sparse,
      row-oriented matrices.}
  }
}
\examples{
## A random sparse matrix :
set.seed(7)
m <- matrix(0, 5, 5)
m[sample(length(m), size = 14)] <- rep(1:9, length=14)
(mm <- as(m, "dgCMatrix"))

tril(mm)        # lower triangle
tril(mm, -1)    # strict lower triangle
triu(mm,  1)    # strict upper triangle
band(mm, -1, 2) # general band
\dontshow{ ## this uses special methods
(x.x <- crossprod(mm))
tril(x.x)
xx <- tril(x.x) + triu(x.x, 1) ## the same as x.x (but stored differently):
txx <- t(as(xx, "dsCMatrix"))
stopifnot(identical(triu(x.x), t(tril(x.x))),
	  identical(class(x.x), class(txx)),
	  identical(as(x.x, "dgCMatrix"), as(txx, "dgCMatrix")))
}
}
\keyword{methods}
\keyword{algebra}
