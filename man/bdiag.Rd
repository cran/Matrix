\name{bdiag}
\alias{bdiag}
\alias{.bdiag}
\title{Construct a Block Diagonal Matrix}
\description{
  Build a block diagonal matrix given several building block matrices.
}
\usage{
bdiag(\dots)
.bdiag(lst)
}
\arguments{
  \item{\dots}{individual matrices or a \code{\link{list}} of matrices.}
  \item{lst}{non-empty \code{\link{list}} of matrices.}
}
\details{
  For non-trivial argument list, \code{bdiag()} calls \code{.bdiag()}.
  The latter maybe useful to programmers.
}
\note{This function has been written and is efficient for the case of
  relatively few block matrices which are typically sparse themselves.

  It is currently \emph{inefficient} for the case of many small dense
  block matrices.
  For the case of \emph{many} dense \eqn{k \times k}{k * k} matrices,
  the \code{bdiag_m()} function in the \sQuote{Examples} is an order of
  magnitude faster.
}
\value{
  A \emph{sparse} matrix obtained by combining the arguments into a
  block diagonal matrix.

  The value of \code{bdiag()} inherits from class
  \code{\linkS4class{CsparseMatrix}}, whereas
  \code{.bdiag()} returns a \code{\linkS4class{TsparseMatrix}}.
}
\author{Martin Maechler, built on a version posted by Berton Gunter to
  R-help; earlier versions have been posted by other authors, notably
  Scott Chasalow to S-news.  Doug Bates's faster implementation builds
  on \code{\linkS4class{TsparseMatrix}} objects.
}
\seealso{\code{\link{Diagonal}} for constructing matrices of
  class \code{\linkS4class{diagonalMatrix}}, or \code{\link{kronecker}}
  which also works for \code{"Matrix"} inheriting matrices.

  \code{\link{bandSparse}} constructs a \emph{banded} sparse matrix from
  its non-zero sub-/super - diagonals.

  Note that other CRAN \R packages have own versions of \code{bdiag()}
  which return traditional matrices.
}
\examples{
bdiag(matrix(1:4, 2), diag(3))
## combine "Matrix" class and traditional matrices:
bdiag(Diagonal(2), matrix(1:3, 3,4), diag(3:2))

mlist <- list(1, 2:3, diag(x=5:3), 27, cbind(1,3:6), 100:101)
bdiag(mlist)
stopifnot(identical(bdiag(mlist), % <- used to fail in earlier versions
                    bdiag(lapply(mlist, as.matrix))))

ml <- c(as(matrix((1:24)\%\% 11 == 0, 6,4),"nMatrix"),
        rep(list(Diagonal(2, x=TRUE)), 3))
mln <- c(ml, Diagonal(x = 1:3))
stopifnot(is(bdiag(ml), "lsparseMatrix"),% failed in Matrix <= 1.0-2
          is(bdiag(mln),"dsparseMatrix") )

## random (diagonal-)block-triangular matrices:
rblockTri <- function(nb, max.ni, lambda = 3) {
   .bdiag(replicate(nb, {
         n <- sample.int(max.ni, 1)
         tril(Matrix(rpois(n*n, lambda=lambda), n,n)) }))
}

(T4 <- rblockTri(4, 10, lambda = 1))
image(T1 <- rblockTri(12, 20))


##' Fast version of Matrix :: .bdiag() -- for the case of *many*  (k x k) matrices:
##' @param lmat list(<mat1>, <mat2>, ....., <mat_N>)  where each mat_j is a  k x k 'matrix'
##' @return a sparse (N*k x N*k) matrix of class  \code{"\linkS4class{dgCMatrix}"}.
bdiag_m <- function(lmat) {
    ## Copyright (C) 2016 Martin Maechler, ETH Zurich
    if(!length(lmat)) return(new("dgCMatrix"))
    stopifnot(is.list(lmat), is.matrix(lmat[[1]]),
              (k <- (d <- dim(lmat[[1]]))[1]) == d[2], # k x k
              all(vapply(lmat, dim, integer(2)) == k)) # all of them
    N <- length(lmat)
    if(N * k > .Machine$integer.max)
        stop("resulting matrix too large; would be  M x M, with M=", N*k)
    M <- as.integer(N * k)
    ## result: an   M x M  matrix
    new("dgCMatrix", Dim = c(M,M),
        ## 'i :' maybe there's a faster way (w/o matrix indexing), but elegant?
        i = as.vector(matrix(0L:(M-1L), nrow=k)[, rep(seq_len(N), each=k)]),
        p = k * 0L:M,
        x = as.double(unlist(lmat, recursive=FALSE, use.names=FALSE)))
}

l12 <- replicate(12, matrix(rpois(16, lambda = 6.4), 4,4), simplify=FALSE)
dim(T12 <- bdiag_m(l12))# 48 x 48
T12[1:20, 1:20]
}
\keyword{array}
