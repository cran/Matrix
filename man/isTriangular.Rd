\name{isTriangular-methods}
\title{Test whether a Matrix is Triangular or Diagonal}
%
\docType{methods}
\keyword{array}
\keyword{programming}
\keyword{methods}
%
\alias{isTriangular}
\alias{isTriangular-methods}
\alias{isDiagonal}
\alias{isDiagonal-methods}
%
\alias{isTriangular,diagonalMatrix-method}
\alias{isTriangular,indMatrix-method}
\alias{isTriangular,matrix-method}
\alias{isTriangular,symmetricMatrix-method}
\alias{isTriangular,triangularMatrix-method}
%
\alias{isDiagonal,CsparseMatrix-method}
\alias{isDiagonal,RsparseMatrix-method}
\alias{isDiagonal,TsparseMatrix-method}
\alias{isDiagonal,diagonalMatrix-method}
\alias{isDiagonal,indMatrix-method}
\alias{isDiagonal,matrix-method}
\alias{isDiagonal,packedMatrix-method}
\alias{isDiagonal,unpackedMatrix-method}
% disambiguation:
\alias{isTriangular,dgCMatrix-method}
\alias{isTriangular,dgRMatrix-method}
\alias{isTriangular,dgTMatrix-method}
\alias{isTriangular,dgeMatrix-method}
\alias{isTriangular,lgCMatrix-method}
\alias{isTriangular,lgRMatrix-method}
\alias{isTriangular,lgTMatrix-method}
\alias{isTriangular,lgeMatrix-method}
\alias{isTriangular,ngCMatrix-method}
\alias{isTriangular,ngRMatrix-method}
\alias{isTriangular,ngTMatrix-method}
\alias{isTriangular,ngeMatrix-method}
%
\description{
  \code{isTriangular} and \code{isDiagonal} test whether their argument
  is a triangular or diagonal matrix, respectively.  Unlike the analogous
  \code{\link{isSymmetric}}, these two functions are generically
  from \pkg{Matrix} rather than \code{base}.  Hence \pkg{Matrix}
  defines methods for traditional matrices of implicit \link{class}
  \code{"\link{matrix}"} in addition to matrices inheriting from
  virtual class \code{"\linkS4class{Matrix}"}.

  By our definition, triangular and diagonal matrices are \emph{square},
  i.e., they have the same number of rows and columns.
}
\usage{
isTriangular(object, upper = NA, \dots)

isDiagonal(object)
}
\arguments{
  \item{object}{an \R object, typically a matrix.}
  \item{upper}{a \link{logical}, either \code{TRUE} or \code{FALSE},
    in which case \code{TRUE} is returned only for upper or lower
    triangular \code{object}; or otherwise \code{\link{NA}} (the
    default), in which case \code{TRUE} is returned for any triangular
    \code{object}.}
  \item{\dots}{further arguments passed to methods
    (currently unused by \pkg{Matrix}).}
}
\value{
  A \link{logical}, either \code{TRUE} or \code{FALSE}
  (never \code{\link{NA}}).

  If \code{object} is triangular and \code{upper} is \code{NA}, then
  \code{isTriangular} returns \code{TRUE} with an \link[=attr]{attribute} 
  \code{kind}, either \code{"U"} or \code{"L"}, indicating that
  \code{object} is \bold{u}pper or \bold{l}ower triangular, respectively.  
  Users should not rely on how \code{kind} is determined for diagonal
  matrices, which are both upper and lower triangular.
}
\seealso{
  \code{\link{isSymmetric}};
  virtual classes \code{"\linkS4class{triangularMatrix}"} and
  \code{"\linkS4class{diagonalMatrix}"} and their subclasses.
}
\examples{
isTriangular(Diagonal(4))
## is TRUE: a diagonal matrix is also (both upper and lower) triangular
(M <- Matrix(c(1,2,0,1), 2,2))
isTriangular(M) # TRUE (*and* of formal class "dtrMatrix")
isTriangular(as(M, "generalMatrix")) # still triangular, even if not "formally"
isTriangular(crossprod(M)) # FALSE

isDiagonal(matrix(c(2,0,0,1), 2,2)) # TRUE

## Look at implementations:
showMethods("isTriangular", includeDefs = TRUE)
showMethods("isDiagonal", includeDefs = TRUE)
}
