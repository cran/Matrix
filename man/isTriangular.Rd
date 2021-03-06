\name{isTriangular}
\title{isTriangular() and isDiagonal() Methods}
%\docType{methods} both generic *and* methods
\alias{isDiagonal}
\alias{isDiagonal-methods}
\alias{isDiagonal,symmetricMatrix-method}
\alias{isDiagonal,triangularMatrix-method}
\alias{isDiagonal,denseMatrix-method}
\alias{isDiagonal,diagonalMatrix-method}
\alias{isDiagonal,sparseMatrix-method}
\alias{isDiagonal,CsparseMatrix-method}
\alias{isDiagonal,TsparseMatrix-method}
\alias{isDiagonal,matrix-method}
%
\alias{isTriangular}
\alias{isTriangular-methods}
\alias{isTriangular,triangularMatrix-method}
\alias{isTriangular,denseMatrix-method}
\alias{isTriangular,diagonalMatrix-method}
\alias{isTriangular,BunchKaufman-method}
\alias{isTriangular,Cholesky-method}
\alias{isTriangular,CsparseMatrix-method}
\alias{isTriangular,TsparseMatrix-method}
\alias{isTriangular,dtCMatrix-method}
\alias{isTriangular,dtRMatrix-method}
\alias{isTriangular,dtTMatrix-method}
\alias{isTriangular,dtpMatrix-method}
\alias{isTriangular,dtrMatrix-method}
\alias{isTriangular,ltCMatrix-method}
\alias{isTriangular,ltRMatrix-method}
\alias{isTriangular,ltTMatrix-method}
\alias{isTriangular,ltpMatrix-method}
\alias{isTriangular,ltrMatrix-method}
\alias{isTriangular,ntCMatrix-method}
\alias{isTriangular,ntRMatrix-method}
\alias{isTriangular,ntTMatrix-method}
\alias{isTriangular,ntpMatrix-method}
\alias{isTriangular,ntrMatrix-method}
\alias{isTriangular,pBunchKaufman-method}
\alias{isTriangular,pCholesky-method}
\alias{isTriangular,matrix-method}
%
\description{
  \code{isTriangular(M)} returns a \code{\link{logical}} indicating
  if \code{M} is a triangular matrix.  Analogously,
  \code{isDiagonal(M)} is true iff \code{M} is a diagonal matrix.

  Contrary to \code{\link{isSymmetric}()}, these two functions are
  generically from package \pkg{Matrix}, and hence also define methods
  for traditional (\code{\link{class}} \code{"matrix"}) matrices.

  By our definition, triangular, diagonal and symmetric matrices are all
  \emph{square}, i.e. have the same number of rows and columns.
}
\usage{
isDiagonal(object)

isTriangular(object, upper = NA, \dots)
}
\arguments{
  \item{object}{any \R object, typically a matrix (traditional or Matrix
    package).}
  \item{upper}{logical, one of \code{NA} (default), \code{FALSE}, or
    \code{TRUE} where the last two cases require a lower or upper
    triangular \code{object} to result in \code{TRUE}.}
  \item{\dots}{potentially further arguments for other methods.}
}
\value{
  a (\dQuote{scalar}) logical, \code{TRUE} or \code{FALSE}, never
  \code{\link{NA}}.  For \code{isTriangular()}, if the result is
  \code{TRUE}, it may contain an attribute (see \code{\link{attributes}}
  \code{"kind"}, either   \code{"L"} or \code{"U"} indicating if it is a
  \bold{l}ower or \bold{u}pper triangular matrix.
}
\seealso{
  \code{\link{isSymmetric}}; formal class (and subclasses)
  \code{"\linkS4class{triangularMatrix}"} and
  \code{"\linkS4class{diagonalMatrix}"}.
}
\examples{
isTriangular(Diagonal(4))
## is TRUE: a diagonal matrix is also (both upper and lower) triangular
(M <- Matrix(c(1,2,0,1), 2,2))
isTriangular(M) # TRUE (*and* of formal class "dtrMatrix")
isTriangular(as(M, "dgeMatrix")) # still triangular, even if not "formally"
isTriangular(crossprod(M)) # FALSE

isDiagonal(matrix(c(2,0,0,1), 2,2)) # TRUE
}
\keyword{methods}
