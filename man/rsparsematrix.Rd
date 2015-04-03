\name{rsparsematrix}
\alias{rsparsematrix}
\title{Random Sparse Matrix}
\description{
  Generate a Random Sparse Matrix Efficiently.
}
\usage{
rsparsematrix(nrow, ncol, density, nnz = round(density * maxE),
              symmetric = FALSE,
              rand.x = function(n) signif(rnorm(nnz), 2), \dots)
}
\arguments{
  \item{nrow, ncol}{number of rows and columns, i.e., the matrix
    dimension (\code{\link{dim}}).}
  \item{density}{optional number in \eqn{[0,1]}, the density is the
    proportion of non-zero entries among all matrix entries.  If
    specified it determines the default for \code{nnz}, otherwise
    \code{nnz} needs to be specified.}
  \item{nnz}{number of non-zero entries, for a sparse matrix typically
    considerably smaller than \code{nrow*ncol}.  Must be specified if
    \code{density} is not.}
  \item{symmetric}{logical indicating if result should be a matrix of
    class \code{\linkS4class{symmetricMatrix}}.  Note that in the symmetric
    case, \code{nnz} denotes the number of non zero entries of the upper
    (or lower) part of the matrix, including the diagonal.}
  \item{rand.x}{\code{\link{NULL}} or the random number generator for the \code{x} slot, a
    \code{\link{function}} such that \code{rand.x(n)} generates a
    numeric vector of length \code{n}.  Typical examples are
    \code{rand.x = rnorm}, or \code{rand.x = runif}; the default is nice
    for didactical purposes.}
  \item{\dots}{optionally further arguments passed to
    \code{\link{sparseMatrix}()}, notably \code{giveCsparse}.}
}
\details{
  The algorithm first samples \dQuote{encoded} \eqn{(i,j)}s without
  replacement, via one dimensional indices, if not \code{symmetric}
  \code{\link{sample.int}(nrow*ncol, nnz)}, then---if \code{rand.x} is
  not \code{NULL}---gets \code{x <- rand.x(nnz)} and calls
  \code{\link{sparseMatrix}(i=i, j=j, x=x, ..)}.  When
  \code{rand.x=NULL}, \code{\link{sparseMatrix}(i=i, j=j, ..)} will
  return a patter\bold{n} matrix (i.e., inheriting from
  \code{\linkS4class{nsparseMatrix}}).
}
\value{
  a \code{\linkS4class{sparseMatrix}}, say \code{M} of dimension (nrow,
  ncol), i.e., with \code{dim(M) == c(nrow, ncol)}, if \code{symmetric}
  is not true, with \code{nzM <- \link{nnzero}(M)} fulfilling
  \code{nzM <= nnz} and typically, \code{nzM == nnz}.
}
\author{Martin Maechler}
\examples{
set.seed(17)# to be reproducible
M <- rsparsematrix(8, 12, nnz = 30) # small example, not very sparse
M
M1 <- rsparsematrix(1000, 20,  nnz = 123,  rand.x = runif)
summary(M1)

## a random *symmetric* Matrix
(S9 <- rsparsematrix(9, 9, nnz = 10, symmetric=TRUE)) # dsCMatrix
nnzero(S9)# ~ 20: as 'nnz' only counts one "triangle"

## a random patter*n* aka boolean Matrix (no 'x' slot):
(n7 <- rsparsematrix(5, 12, nnz = 10, rand.x = NULL))

## a [T]riplet representation sparseMatrix:
T2 <- rsparsematrix(40, 12, nnz = 99, giveCsparse=FALSE)
head(T2)
}
\keyword{array}
\keyword{distribution}

