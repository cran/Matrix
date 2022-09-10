\name{[-methods}
\title{Methods for "[": Extraction or Subsetting in Package 'Matrix'}
\docType{methods}
%
\alias{[-methods}
%
\alias{[,CsparseMatrix,index,index,logical-method}
\alias{[,CsparseMatrix,index,missing,logical-method}
\alias{[,CsparseMatrix,missing,index,logical-method}
\alias{[,Matrix,ANY,ANY,ANY-method}
\alias{[,Matrix,index,index,missing-method}
\alias{[,Matrix,index,missing,missing-method}
\alias{[,Matrix,lMatrix,missing,missing-method}
\alias{[,Matrix,logical,missing,missing-method}
\alias{[,Matrix,matrix,missing,ANY-method}
\alias{[,Matrix,matrix,missing,missing-method}
\alias{[,Matrix,missing,index,missing-method}
\alias{[,Matrix,missing,missing,ANY-method}
\alias{[,Matrix,missing,missing,logical-method}
\alias{[,Matrix,missing,missing,missing-method}
\alias{[,Matrix,nMatrix,missing,missing-method}
\alias{[,TsparseMatrix,index,index,logical-method}
\alias{[,TsparseMatrix,index,missing,logical-method}
\alias{[,TsparseMatrix,missing,index,logical-method}
\alias{[,abIndex,index,ANY,ANY-method}
\alias{[,denseMatrix,index,index,logical-method}
\alias{[,denseMatrix,index,missing,logical-method}
\alias{[,denseMatrix,matrix,missing,ANY-method}
\alias{[,denseMatrix,matrix,missing,missing-method}
\alias{[,denseMatrix,missing,index,logical-method}
\alias{[,diagonalMatrix,index,index,logical-method}
\alias{[,diagonalMatrix,index,missing,logical-method}
\alias{[,diagonalMatrix,missing,index,logical-method}
\alias{[,indMatrix,index,missing,logical-method}
\alias{[,packedMatrix,NULL,NULL,logical-method}
\alias{[,packedMatrix,NULL,NULL,missing-method}
\alias{[,packedMatrix,NULL,index,logical-method}
\alias{[,packedMatrix,NULL,index,missing-method}
\alias{[,packedMatrix,NULL,missing,logical-method}
\alias{[,packedMatrix,NULL,missing,missing-method}
\alias{[,packedMatrix,index,NULL,logical-method}
\alias{[,packedMatrix,index,NULL,missing-method}
\alias{[,packedMatrix,index,index,logical-method}
\alias{[,packedMatrix,index,index,missing-method}
\alias{[,packedMatrix,index,missing,logical-method}
\alias{[,packedMatrix,index,missing,missing-method}
\alias{[,packedMatrix,lMatrix,NULL,logical-method}
\alias{[,packedMatrix,lMatrix,NULL,missing-method}
\alias{[,packedMatrix,lMatrix,index,logical-method}
\alias{[,packedMatrix,lMatrix,index,missing-method}
\alias{[,packedMatrix,lMatrix,missing,logical-method}
\alias{[,packedMatrix,lMatrix,missing,missing-method}
\alias{[,packedMatrix,matrix,NULL,logical-method}
\alias{[,packedMatrix,matrix,NULL,missing-method}
\alias{[,packedMatrix,matrix,index,logical-method}
\alias{[,packedMatrix,matrix,index,missing-method}
\alias{[,packedMatrix,matrix,missing,logical-method}
\alias{[,packedMatrix,matrix,missing,missing-method}
\alias{[,packedMatrix,missing,NULL,logical-method}
\alias{[,packedMatrix,missing,NULL,missing-method}
\alias{[,packedMatrix,missing,index,logical-method}
\alias{[,packedMatrix,missing,index,missing-method}
\alias{[,packedMatrix,missing,missing,logical-method}
\alias{[,packedMatrix,missing,missing,missing-method}
\alias{[,sparseMatrix,index,index,logical-method}
\alias{[,sparseMatrix,index,missing,logical-method}
\alias{[,sparseMatrix,missing,index,logical-method}
\alias{[,sparseVector,index,ANY,ANY-method}
\alias{[,sparseVector,lsparseVector,ANY,ANY-method}
\alias{[,sparseVector,nsparseVector,ANY,ANY-method}
%
\description{
  Methods for \code{"["}, i.e., extraction or subsetting mostly of
  matrices, in package \pkg{Matrix}.
}
\section{Methods}{
  There are more than these:
  \describe{
    \item{x = "Matrix", i = "missing", j = "missing", drop= "ANY"}{ ... }
    \item{x = "Matrix", i = "numeric", j = "missing", drop= "missing"}{ ... }
    \item{x = "Matrix", i = "missing", j = "numeric", drop= "missing"}{ ... }

    \item{x = "dsparseMatrix", i = "missing", j = "numeric", drop= "logical"}{ ... }
    \item{x = "dsparseMatrix", i = "numeric", j = "missing", drop= "logical"}{ ... }
    \item{x = "dsparseMatrix", i = "numeric", j = "numeric", drop= "logical"}{ ... }
  }
}
\seealso{
  \code{\link{[<--methods}} for sub\emph{assign}ment to \code{"Matrix"}
  objects.
  \code{\link{Extract}} about the standard extraction.
}
\examples{
% regression tests are in ../tests/indexing.R
str(m <- Matrix(round(rnorm(7*4),2), nrow = 7))
stopifnot(identical(m, m[]))
m[2, 3]   # simple number
m[2, 3:4] # simple numeric of length 2
m[2, 3:4, drop=FALSE] # sub matrix of class 'dgeMatrix'
## rows or columns only:
m[1,]     # first row, as simple numeric vector
m[,1:2]   # sub matrix of first two columns

showMethods("[", inherited = FALSE)
}
\keyword{methods}
\keyword{array}
