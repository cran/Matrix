\name{cBind}
\title{'cbind()' and 'rbind()' recursively built on cbind2/rbind2}
\alias{cBind}
\alias{rBind}
%no \alias{cbind2} or \alias{rbind2}
% many more cbind2|rbind2 methods are in *-class.Rd, e.g. ./sparseMatrix-class.Rd
\alias{cbind2,denseMatrix,sparseMatrix-method}
\alias{cbind2,sparseMatrix,denseMatrix-method}
\alias{rbind2,denseMatrix,sparseMatrix-method}
\alias{rbind2,sparseMatrix,denseMatrix-method}
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
\section{Historical Remark}{
  Before \R version 3.2.0 (April 2015), we have needed a substitute for
  \emph{S4-enabled} versions of \code{cbind} and \code{rbind}, and
  provided \code{cBind} and \code{rBind} with identical syntax and
  semantic in order to bind together multiple matrices (\code{"matrix"}
  or \code{"Matrix"} and vectors.
  With \R version 3.2.0 and newer, \code{cBind} and \code{rBind} are
  \emph{deprecated} and produce a deprecation warning (via
  \code{\link{.Deprecated}}), and your code should start using
  \code{cbind()} and \code{rbind()} instead.
}
\usage{
## cbind(..., deparse.level = 1)
## rbind(..., deparse.level = 1)

## and e.g.,
\S4method{cbind2}{denseMatrix,sparseMatrix}(x,y, sparse = NA, ...)
\S4method{cbind2}{sparseMatrix,denseMatrix}(x,y, sparse = NA, ...)
\S4method{rbind2}{denseMatrix,sparseMatrix}(x,y, sparse = NA, ...)
\S4method{rbind2}{sparseMatrix,denseMatrix}(x,y, sparse = NA, ...)
}
\arguments{
  \item{\dots, x, y}{matrix-like \R objects to be bound together, see
    \code{\link{cbind}} and \code{\link{rbind}}.}
  \item{sparse}{option \code{\link{logical}} indicating if the result
    should be sparse, i.e., formally inheriting from \code{"\linkS4class{sparseMatrix}"}.
    The default, \code{\link{NA}}, decides from the \dQuote{sparsity} of
    \code{x} and \code{y}, see e.g., the \R code in
    \code{selectMethod(cbind2, c("sparseMatrix","denseMatrix"))}.}
  \item{deparse.level}{integer determining under which circumstances
    column and row names are built from the actual arguments'
    \sQuote{expression}, see \code{\link{cbind}}.}
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
\seealso{\code{\link{cbind2}}, \code{\link{cbind}},
  Documentation in base \R's \pkg{methods} package.

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
\keyword{array}
\keyword{manip}

