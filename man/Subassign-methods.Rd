\name{[<--methods}
\docType{methods}
\alias{[<--methods}
\alias{Subassign-methods}% <- nicer to use in \link{.}
%
\title{Methods for "[<-" - Assigning to Subsets for 'Matrix'}
% in ../R/Matrix.R
\alias{[<-,Matrix,ANY,ANY,ANY-method}
\alias{[<-,Matrix,ANY,ANY,Matrix-method}
\alias{[<-,Matrix,ANY,ANY,matrix-method}
\alias{[<-,Matrix,missing,missing,ANY-method}
% \alias{[<-,Matrix,missing,numeric,missing-method}
% \alias{[<-,Matrix,numeric,missing,missing-method}
% \alias{[<-,Matrix,numeric,numeric,missing-method}
% in ../R/denseMatrix.R
\alias{[<-,denseMatrix,index,missing,numeric-method}
\alias{[<-,denseMatrix,index,index,numeric-method}
\alias{[<-,denseMatrix,missing,index,numeric-method}
% in ../R/diagMatrix.R
\alias{[<-,diagonalMatrix,ANY,ANY,ANY-method}
% -> ../R/Tsparse.R : replValue = { numeric, logical }
\alias{[<-,TsparseMatrix,index,index,replValue-method}
\alias{[<-,TsparseMatrix,index,missing,replValue-method}
\alias{[<-,TsparseMatrix,missing,index,replValue-method}
% \alias{[<-,dgTMatrix,missing,missing,numeric-method}
\alias{[<-,CsparseMatrix,index,index,replValue-method}
\alias{[<-,CsparseMatrix,index,missing,replValue-method}
\alias{[<-,CsparseMatrix,missing,index,replValue-method}
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
  %% ./Xtrct-methods.Rd:
  \code{\link{[-methods}} for subsetting \code{"Matrix"} objects; the
  \code{\linkS4class{index}} class;
  \code{\link{Extract}} about the standard subset assignment (and extraction).
}
\examples{
%% Note that ./Xtrct-methods.Rd has the indexing ones
set.seed(101)
(a <- m <- Matrix(round(rnorm(7*4),2), nrow = 7))
% not yet available%% FIXME!
% a[] <- 2.2
m[2, 3]  <- 3.14 # simple number
m[3, 3:4]<- 3:4  # simple numeric of length 2

## sub matrix assignment:
m[-(4:7), 3:4] <- cbind(1,2:4) #-> upper right corner of 'm'
m[3:5, 2:3] <- 0
m[6:7, 1:2] <- Diagonal(2)
m

## rows or columns only:
m[1,] <- 10
m[,2] <- 1:7
m[-(1:6), ] <- 3:0 # not the first 6 rows, i.e. only the 7th
as(m, "sparseMatrix")
}
\keyword{methods}
\keyword{array}
