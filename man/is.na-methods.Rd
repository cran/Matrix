\name{is.na-methods}
\title{is.na(), is.finite() Methods for 'Matrix' Objects}
%
\docType{methods}
\keyword{NA}
\keyword{math}
\keyword{programming}
\keyword{methods}
%
\alias{is.na}
\alias{is.na-methods}
\alias{is.nan}
\alias{is.nan-methods}
\alias{is.finite}
\alias{is.finite-methods}
\alias{is.infinite}
\alias{is.infinite-methods}
\alias{anyNA}
\alias{anyNA-methods}
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
\alias{is.nan,ddiMatrix-method}
\alias{is.nan,dgeMatrix-method}
\alias{is.nan,dspMatrix-method}
\alias{is.nan,dsparseMatrix-method}
\alias{is.nan,dsyMatrix-method}
\alias{is.nan,dtpMatrix-method}
\alias{is.nan,dtrMatrix-method}
\alias{is.nan,indMatrix-method}
\alias{is.nan,lMatrix-method}
\alias{is.nan,nMatrix-method}
\alias{is.nan,nsparseVector-method}
\alias{is.nan,sparseVector-method}
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
\description{
  Methods for generic functions \code{\link{is.na}()},
  \code{\link{is.nan}()}, \code{\link{is.finite}()},
  \code{\link{is.infinite}()}, and \code{\link{anyNA}()},
  for objects inheriting from virtual class
  \code{\linkS4class{Matrix}} or \code{\linkS4class{sparseVector}}.
}
\usage{
\S4method{is.na}{dsparseMatrix}(x)
\S4method{is.nan}{dsparseMatrix}(x)
\S4method{is.finite}{dsparseMatrix}(x)
\S4method{is.infinite}{dsparseMatrix}(x)
\S4method{anyNA}{dsparseMatrix}(x)
## ...
## and for other classes
}
\arguments{
  \item{x}{an \R object, here a sparse or dense matrix or vector.}
}
\value{
  For \code{is.*()}, an \code{\linkS4class{nMatrix}} or
  \code{\linkS4class{nsparseVector}} matching the dimensions
  of \code{x} and specifying the positions in \code{x} of
  (some subset of) \code{\link{NA}}, \code{\link{NaN}},
  \code{\link{Inf}}, and \code{-Inf}.
  For \code{anyNA}(), \code{TRUE} if \code{x} contains \code{NA}
  or \code{NaN} and \code{FALSE} otherwise.
}
\seealso{
  \code{\link{NA}}, \code{\link{NaN}}, \code{\link{Inf}}
}
\examples{
(M <- Matrix(1:6, nrow = 4, ncol = 3,
             dimnames = list(letters[1:4], LETTERS[1:3])))
stopifnot(!anyNA(M), !any(is.na(M)))

M[2:3, 2] <- NA
(inM <- is.na(M))
stopifnot(anyNA(M), sum(inM) == 2)

(A <- spMatrix(nrow = 10, ncol = 20,
               i = c(1, 3:8), j = c(2, 9, 6:10), x = 7 * (1:7)))
stopifnot(!anyNA(A), !any(is.na(A)))

A[2, 3] <- A[1, 2] <- A[5, 5:9] <- NA
(inA <- is.na(A))
stopifnot(anyNA(A), sum(inA) == 1 + 1 + 5)
}
