\name{is.na-methods}
\title{is.na(), is.{in}finite() Methods for 'Matrix' Objects}
\docType{methods}
%
\alias{is.na-methods}
%
\alias{anyNA,ddenseMatrix-method}
\alias{anyNA,diagonalMatrix-method}
\alias{anyNA,dsparseMatrix-method}
\alias{anyNA,indMatrix-method}
\alias{anyNA,ldenseMatrix-method}
\alias{anyNA,lsparseMatrix-method}
\alias{anyNA,nMatrix-method}
\alias{anyNA,nsparseVector-method}
\alias{anyNA,sparseVector-method}
%
\alias{is.na,abIndex-method}
\alias{is.na,dgeMatrix-method}
\alias{is.na,diagonalMatrix-method}
\alias{is.na,dspMatrix-method}
\alias{is.na,dsparseMatrix-method}
\alias{is.na,dsyMatrix-method}
\alias{is.na,dtpMatrix-method}
\alias{is.na,dtrMatrix-method}
\alias{is.na,indMatrix-method}
\alias{is.na,lgeMatrix-method}
\alias{is.na,lspMatrix-method}
\alias{is.na,lsparseMatrix-method}
\alias{is.na,lsyMatrix-method}
\alias{is.na,ltpMatrix-method}
\alias{is.na,ltrMatrix-method}
\alias{is.na,nMatrix-method}
\alias{is.na,nsparseVector-method}
\alias{is.na,sparseVector-method}
%
\alias{is.finite,abIndex-method}
\alias{is.finite,dgeMatrix-method}
\alias{is.finite,diagonalMatrix-method}
\alias{is.finite,dspMatrix-method}
\alias{is.finite,dsparseMatrix-method}
\alias{is.finite,dsyMatrix-method}
\alias{is.finite,dtpMatrix-method}
\alias{is.finite,dtrMatrix-method}
\alias{is.finite,indMatrix-method}
\alias{is.finite,lgeMatrix-method}
\alias{is.finite,lspMatrix-method}
\alias{is.finite,lsparseMatrix-method}
\alias{is.finite,lsyMatrix-method}
\alias{is.finite,ltpMatrix-method}
\alias{is.finite,ltrMatrix-method}
\alias{is.finite,nMatrix-method}
\alias{is.finite,nsparseVector-method}
\alias{is.finite,sparseVector-method}
%
\alias{is.infinite,abIndex-method}
\alias{is.infinite,ddiMatrix-method}
\alias{is.infinite,dgeMatrix-method}
\alias{is.infinite,dspMatrix-method}
\alias{is.infinite,dsparseMatrix-method}
\alias{is.infinite,dsyMatrix-method}
\alias{is.infinite,dtpMatrix-method}
\alias{is.infinite,dtrMatrix-method}
\alias{is.infinite,indMatrix-method}
\alias{is.infinite,lMatrix-method}
\alias{is.infinite,nMatrix-method}
\alias{is.infinite,nsparseVector-method}
\alias{is.infinite,sparseVector-method}
%
\description{
  Methods for functions \code{\link{is.na}()}, \code{\link{is.finite}()},
  and \code{\link{is.infinite}()} for all Matrix objects
  (those extending class \code{\linkS4class{Matrix}}):

  \describe{
    \item{x = "denseMatrix"}{returns an \code{\linkS4class{nMatrix}}
      of the same dimension as \code{x}, with TRUE whenever \code{x}
      is \code{\link{NA}}, finite, or infinite, respectively.}
    \item{x = "sparseMatrix"}{ditto.}
  }
}
\usage{
\S4method{is.na}{lsparseMatrix}(x)
\S4method{is.finite}{dsparseMatrix}(x)
\S4method{is.infinite}{dgeMatrix}(x)
## ...
## and for other classes

\S4method{anyNA}{nMatrix}(x)
## ...
## and for other classes
}
\arguments{
  \item{x}{sparse or dense matrix or sparse vector (here; any \R object
    in general).}
}
\seealso{
  \code{\link{NA}}, \code{\link{is.na}};
  \code{\link{is.finite}}, \code{\link{is.infinite}};
  \code{\linkS4class{nMatrix}},
  \code{\linkS4class{denseMatrix}}, \code{\linkS4class{sparseMatrix}}.

  The \code{\linkS4class{sparseVector}} class.
}
\examples{
M <- Matrix(1:6, nrow=4, ncol=3,
       dimnames = list(c("a", "b", "c", "d"), c("A", "B", "C")))
stopifnot(all(!is.na(M)))
M[2:3,2] <- NA
is.na(M)
if(exists("anyNA", mode="function"))
anyNA(M)

A <- spMatrix(10,20, i = c(1,3:8),
                     j = c(2,9,6:10),
                     x = 7 * (1:7))
stopifnot(all(!is.na(A)))
%TODO stopifnot(all(!is.na(A)), !anyNA(A))
A[2,3] <- A[1,2] <- A[5, 5:9] <- NA
inA <- is.na(A)
stopifnot(sum(inA) == 1+1+5)
\dontshow{
if(exists("anyNA", mode="function"))
stopifnot(anyNA(M), anyNA(A))
}% dont..
}
\keyword{methods}
