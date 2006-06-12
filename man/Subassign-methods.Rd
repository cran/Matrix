\name{[<--methods}
\docType{methods}
\alias{[<--methods}
\alias{Subassign-methods}% <- nicer to use in \link{.}
%
\title{Methods for "[<-" - Assigning to Subsets for 'Matrix'}
% in ../R/Matrix.R
\alias{[<-,Matrix,ANY,ANY,ANY-method}
\alias{[<-,Matrix,missing,missing,ANY-method}
% \alias{[<-,Matrix,missing,numeric,missing-method}
% \alias{[<-,Matrix,numeric,missing,missing-method}
% \alias{[<-,Matrix,numeric,numeric,missing-method}
% in ../R/denseMatrix.R
\alias{[<-,denseMatrix,index,missing,numeric-method}
\alias{[<-,denseMatrix,index,index,numeric-method}
\alias{[<-,denseMatrix,missing,index,numeric-method}
%
% % -- maybe replace dgT* by Tsparse* ?
\alias{[<-,dgTMatrix,index,index,numeric-method}
\alias{[<-,dgTMatrix,index,missing,numeric-method}
\alias{[<-,dgTMatrix,missing,index,numeric-method}
% \alias{[<-,dgTMatrix,missing,missing,numeric-method}
\alias{[<-,dgCMatrix,index,index,numeric-method}
\alias{[<-,dgCMatrix,index,missing,numeric-method}
\alias{[<-,dgCMatrix,missing,index,numeric-method}
% %
% \alias{[<-,dsparseMatrix,missing,index,numeric-method}
% \alias{[<-,dsparseMatrix,index,missing,numeric-method}
% \alias{[<-,dsparseMatrix,index,index,numeric-method}
% %
% \alias{[<-,lsparseMatrix,missing,numeric,logical-method}
% \alias{[<-,lsparseMatrix,numeric,missing,logical-method}
% \alias{[<-,lsparseMatrix,numeric,numeric,logical-method}
%-------
\description{
  Methods for \code{"[<-"}, i.e., extraction or subsetting mostly of
  matrices, in package \pkg{Matrix}.
}
\section{Methods}{
  There are more than these:
  \describe{
    \item{x = "Matrix", i = "missing", j = "missing", value= "ANY"}{
      is currently a simple fallback method implementation which ensures
      \dQuote{readable} error messages.}
    \item{x = "Matrix", i = "ANY", j = "ANY", value= "ANY"}{ currently
      gives an error }

    \item{x = "denseMatrix", i = "index", j = "missing", value= "numeric"}{ ... }
    \item{x = "denseMatrix", i = "index", j = "index", value= "numeric"}{ ... }
    \item{x = "denseMatrix", i = "missing", j = "index", value= "numeric"}{ ... }
  }
}
\seealso{
  \code{\link{[-methods}} for subsetting \code{"Matrix"} objects; the
  \code{\linkS4class{index}} class;
  \code{\link{Extract}} about the standard subset assignment (and extraction).
}
\examples{
(m <- Matrix(round(rnorm(7*4),2), nrow = 7))
stopifnot(identical(m, m[]))
m[2, 3]   # simple number
m[2, 3:4] # simple numeric of length 2

m[2, 3:4, drop=FALSE] # sub matrix of class 'dgeMatrix'
m[-(4:7), 3:4]        # dito; the upper right corner of 'm'

## rows or columns only:
m[1,]     # first row, as simple numeric vector
m[,2]     # 2nd column
m[,1:2]   # sub matrix of first two columns
m[-(1:6),, drop=FALSE] # not the first 6 rows, i.e. only the 7th
}
\keyword{methods}
\keyword{array}
