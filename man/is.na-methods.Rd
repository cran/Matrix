\name{is.na-methods}
\title{is.na(), is.finite() Methods for 'Matrix' Objects}
%
\docType{methods}
\keyword{NA}
\keyword{math}
\keyword{programming}
\keyword{methods}
%
\alias{anyNA}
\alias{anyNA-methods}
\alias{is.na}
\alias{is.na-methods}
\alias{is.nan}
\alias{is.nan-methods}
\alias{is.infinite}
\alias{is.infinite-methods}
\alias{is.finite}
\alias{is.finite-methods}
%
\alias{anyNA,denseMatrix-method}
\alias{anyNA,diagonalMatrix-method}
\alias{anyNA,indMatrix-method}
\alias{anyNA,sparseMatrix-method}
\alias{anyNA,sparseVector-method}
%
\alias{is.na,abIndex-method}
\alias{is.na,denseMatrix-method}
\alias{is.na,diagonalMatrix-method}
\alias{is.na,indMatrix-method}
\alias{is.na,sparseMatrix-method}
\alias{is.na,sparseVector-method}
%
\alias{is.nan,denseMatrix-method}
\alias{is.nan,diagonalMatrix-method}
\alias{is.nan,indMatrix-method}
\alias{is.nan,sparseMatrix-method}
\alias{is.nan,sparseVector-method}
%
\alias{is.infinite,abIndex-method}
\alias{is.infinite,denseMatrix-method}
\alias{is.infinite,diagonalMatrix-method}
\alias{is.infinite,indMatrix-method}
\alias{is.infinite,sparseMatrix-method}
\alias{is.infinite,sparseVector-method}
%
\alias{is.finite,abIndex-method}
\alias{is.finite,denseMatrix-method}
\alias{is.finite,diagonalMatrix-method}
\alias{is.finite,indMatrix-method}
\alias{is.finite,sparseMatrix-method}
\alias{is.finite,sparseVector-method}
%
\description{
  Methods for generic functions \code{\link{anyNA}()},
  \code{\link{is.na}()}, \code{\link{is.nan}()},
  \code{\link{is.infinite}()}, and \code{\link{is.finite}()},
  for objects inheriting from virtual class
  \code{\linkS4class{Matrix}} or \code{\linkS4class{sparseVector}}.
}
\usage{
\S4method{is.na}{denseMatrix}(x)
\S4method{is.na}{sparseMatrix}(x)
\S4method{is.na}{diagonalMatrix}(x)
\S4method{is.na}{indMatrix}(x)
\S4method{is.na}{sparseVector}(x)
## ...
## and likewise for anyNA, is.nan, is.infinite, is.finite
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
