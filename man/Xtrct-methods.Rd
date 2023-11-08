\name{Subscript-methods}
\title{Methods for "[": Extraction or Subsetting in Package 'Matrix'}
%
\docType{methods}
\keyword{array}
\keyword{methods}
%
\alias{[}
\alias{[-methods}
\alias{Subscript-methods}
%
\alias{[,Matrix,ANY,NULL,ANY-method}
\alias{[,Matrix,NULL,ANY,ANY-method}
\alias{[,Matrix,NULL,NULL,ANY-method}
\alias{[,Matrix,index,index,logical-method}
\alias{[,Matrix,index,index,missing-method}
\alias{[,Matrix,index,missing,logical-method}
\alias{[,Matrix,index,missing,missing-method}
\alias{[,Matrix,lMatrix,missing,missing-method}
\alias{[,Matrix,matrix,missing,missing-method}
\alias{[,Matrix,missing,index,logical-method}
\alias{[,Matrix,missing,index,missing-method}
\alias{[,Matrix,missing,missing,logical-method}
\alias{[,Matrix,missing,missing,missing-method}
\alias{[,Matrix,nMatrix,missing,missing-method}
\alias{[,abIndex,index,ANY,ANY-method}
\alias{[,sparseVector,NULL,ANY,ANY-method}
\alias{[,sparseVector,index,missing,missing-method}
\alias{[,sparseVector,sparseVector,missing,missing-method}
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
\dontshow{ % for R_DEFAULT_PACKAGES=NULL
library(stats, pos = "package:base", verbose = FALSE)
library(utils, pos = "package:base", verbose = FALSE)
}
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
