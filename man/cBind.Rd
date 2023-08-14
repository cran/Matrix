\name{cbind2-methods}
\title{'cbind()' and 'rbind()' recursively built on cbind2/rbind2}
%
\docType{methods}
\keyword{array}
\keyword{manip}
\keyword{methods}
%
\alias{cbind2}
\alias{cbind2-methods}
\alias{rbind2}
\alias{rbind2-methods}
%
\alias{cbind2,Matrix,Matrix-method}
\alias{cbind2,Matrix,NULL-method}
\alias{cbind2,Matrix,missing-method}
\alias{cbind2,Matrix,vector-method}
\alias{cbind2,NULL,Matrix-method}
\alias{cbind2,ddiMatrix,matrix-method}
\alias{cbind2,ddiMatrix,vector-method}
\alias{cbind2,denseMatrix,denseMatrix-method}
\alias{cbind2,denseMatrix,matrix-method}
\alias{cbind2,denseMatrix,numeric-method}
\alias{cbind2,denseMatrix,sparseMatrix-method}
\alias{cbind2,diagonalMatrix,sparseMatrix-method}
\alias{cbind2,indMatrix,indMatrix-method}
\alias{cbind2,ldiMatrix,matrix-method}
\alias{cbind2,ldiMatrix,vector-method}
\alias{cbind2,matrix,ddiMatrix-method}
\alias{cbind2,matrix,denseMatrix-method}
\alias{cbind2,matrix,ldiMatrix-method}
\alias{cbind2,matrix,sparseMatrix-method}
\alias{cbind2,numeric,denseMatrix-method}
\alias{cbind2,sparseMatrix,denseMatrix-method}
\alias{cbind2,sparseMatrix,diagonalMatrix-method}
\alias{cbind2,sparseMatrix,matrix-method}
\alias{cbind2,sparseMatrix,sparseMatrix-method}
\alias{cbind2,vector,Matrix-method}
\alias{cbind2,vector,ddiMatrix-method}
\alias{cbind2,vector,ldiMatrix-method}
%
\alias{rbind2,Matrix,Matrix-method}
\alias{rbind2,Matrix,NULL-method}
\alias{rbind2,Matrix,missing-method}
\alias{rbind2,Matrix,vector-method}
\alias{rbind2,NULL,Matrix-method}
\alias{rbind2,ddiMatrix,matrix-method}
\alias{rbind2,ddiMatrix,vector-method}
\alias{rbind2,denseMatrix,denseMatrix-method}
\alias{rbind2,denseMatrix,matrix-method}
\alias{rbind2,denseMatrix,numeric-method}
\alias{rbind2,denseMatrix,sparseMatrix-method}
\alias{rbind2,diagonalMatrix,sparseMatrix-method}
\alias{rbind2,indMatrix,indMatrix-method}
\alias{rbind2,ldiMatrix,matrix-method}
\alias{rbind2,ldiMatrix,vector-method}
\alias{rbind2,matrix,ddiMatrix-method}
\alias{rbind2,matrix,denseMatrix-method}
\alias{rbind2,matrix,ldiMatrix-method}
\alias{rbind2,matrix,sparseMatrix-method}
\alias{rbind2,numeric,denseMatrix-method}
\alias{rbind2,sparseMatrix,denseMatrix-method}
\alias{rbind2,sparseMatrix,diagonalMatrix-method}
\alias{rbind2,sparseMatrix,matrix-method}
\alias{rbind2,sparseMatrix,sparseMatrix-method}
\alias{rbind2,vector,Matrix-method}
\alias{rbind2,vector,ddiMatrix-method}
\alias{rbind2,vector,ldiMatrix-method}
%
\description{
  The base functions \code{\link{cbind}} and \code{\link{rbind}} are
  defined for an arbitrary number of arguments and hence have the first
  formal argument \code{...}.  Now, when S4 objects are found among the arguments,
  base \code{cbind()} and \code{rbind()} internally \dQuote{dispatch}
  \emph{recursively}, calling \code{\link{cbind2}} or \code{\link{rbind2}}
  respectively, where these have methods defined and so should dispatch
  appropriately.

  \code{\link{cbind2}()} and \code{\link{rbind2}()} are from the
  \pkg{methods} package, i.e., standard \R, and have been provided for
  binding together \emph{two} matrices, where in \pkg{Matrix}, we have
  defined methods for these and the \code{'Matrix'} matrices.
}
\usage{
## cbind(..., deparse.level = 1)
## rbind(..., deparse.level = 1)

\S4method{cbind2}{denseMatrix,sparseMatrix}(x, y, sparse = NA, \dots)
\S4method{cbind2}{sparseMatrix,denseMatrix}(x, y, sparse = NA, \dots)
\S4method{rbind2}{denseMatrix,sparseMatrix}(x, y, sparse = NA, \dots)
\S4method{rbind2}{sparseMatrix,denseMatrix}(x, y, sparse = NA, \dots)
}
\arguments{
  \item{\dots}{for \code{[cr]bind}, vector- or matrix-like \R objects
	to be bound together; for \code{[cr]bind2}, further arguments
	passed to or from methods; see \code{\link{cbind}} and
	\code{\link[methods]{cbind2}}.}
  \item{deparse.level}{integer controlling the construction of labels
	in the case of non-matrix-like arguments; see \code{\link{cbind}}.}
  \item{x, y}{vector- or matrix-like \R objects to be bound together.}
  \item{sparse}{\code{\link{logical}} indicating if the result should
	be formally sparse, i.e., if it should inherit from virtual class
	\code{\linkS4class{sparseMatrix}}.  \code{\link{NA}}, the default,
	decides based on the \dQuote{sparsity} of \code{x} and \code{y};
	see, e.g., \code{selectMethod(cbind2, c("sparseMatrix", "denseMatrix"))}.}
}
\value{
  typically a \sQuote{matrix-like} object of a similar
  \code{\link{class}} as the first argument in \code{\dots}.

  Note that sometimes by default, the result is a
  \code{\linkS4class{sparseMatrix}} if one of the arguments is (even in
  the case where this is not efficient).  In other cases,
  the result is chosen to be sparse when there are more zero entries is
  than non-zero ones (as the default \code{sparse} in
  \code{\link{Matrix}()}).
}
\author{Martin Maechler}
\seealso{\code{\link{cbind}}, \code{\link[methods]{cbind2}}.

  Our class definition help pages mentioning \code{cbind2()} and
  \code{rbind2()} methods:
  \code{"\linkS4class{denseMatrix}"},
  \code{"\linkS4class{diagonalMatrix}"},
  \code{"\linkS4class{indMatrix}"}.
}
\examples{
(a <- matrix(c(2:1,1:2), 2,2))

(M1 <- cbind(0, rbind(a, 7))) # a traditional matrix

D <- Diagonal(2)
(M2 <- cbind(4, a, D, -1, D, 0)) # a sparse Matrix

stopifnot(validObject(M2), inherits(M2, "sparseMatrix"),
          dim(M2) == c(2,9))
}
