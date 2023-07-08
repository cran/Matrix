\name{dmperm}
\title{Dulmage-Mendelsohn Permutation / Decomposition}
%
\keyword{algebra}
\keyword{array}
\keyword{utilities}
%
\alias{dmperm}
%
\description{
  For any \eqn{n \times m}{n * m} (typically) sparse matrix \code{x}
  compute the Dulmage-Mendelsohn row and columns permutations which at
  first splits the \eqn{n} rows and \code{m} columns into coarse partitions
  each; and then a finer one, reordering rows and columns such that the
  permutated matrix is \dQuote{as upper triangular} as possible.
}
\usage{
dmperm(x, nAns = 6L, seed = 0L)
}
\arguments{
  \item{x}{a typically sparse matrix; internally coerced to either
    \code{"\linkS4class{dgCMatrix}"} or
    \code{"\linkS4class{dtCMatrix}"}.}
  \item{nAns}{an integer specifying the \code{\link{length}} of the
    resulting \code{\link{list}}.  Must be 2, 4, or 6.
    %% FIXME: more
  }
  \item{seed}{an integer code in -1,0,1; determining the (initial)
    permutation; by default, \code{seed = 0}, no (or the identity) permutation;
    \code{seed = -1} uses the \dQuote{reverse} permutation \code{k:1}; for
    \code{seed = 1}, it is a \emph{random} permutation (using R's RNG,
    seed, etc).}
}
\details{
  See the book section by Tim Davis; page 122--127, in the References.
}
\value{
  a named \code{\link{list}} with (by default) 6 components,
 \item{p}{integer vector with the permutation \code{p}, of length \code{nrow(x)}.}
 \item{q}{integer vector with the permutation \code{q}, of length \code{ncol(x)}.}
%% FIXME: give __examples__ below, showing the blocks !!
 \item{r}{integer vector of length \code{nb+1}, where block k is rows r[k] to r[k+1]-1 in A[p,q].}
 \item{s}{integer vector of length \code{nb+1}, where block k is cols s[k] to s[k+1]-1 in A[p,q].}
 \item{rr5}{integer vector of length 5, defining the coarse row
   decomposition.}% FIXME: more specifically ?
 \item{cc5}{integer vector of length 5, defining the coarse column decomposition.}
}
\references{
  Section 7.4 \emph{Dulmage-Mendelsohn decomposition}, pp. 122 ff  of
  \cr
  Timothy A. Davis (2006)
  \emph{Direct Methods for Sparse Linear Systems}, SIAM Series
  \dQuote{Fundamentals of Algorithms}.
}
\author{Martin Maechler, with a lot of \dQuote{encouragement} by Mauricio
  Vargas.
}
%% \note{
%%
%% }
\seealso{
 \code{\link{Schur}}, the class of permutation matrices; \code{"\linkS4class{pMatrix}"}.
}
\examples{
\dontshow{ % for R_DEFAULT_PACKAGES=NULL
library(utils, pos = "package:base", verbose = FALSE)
}
set.seed(17)
(S9 <- rsparsematrix(9, 9, nnz = 10, symmetric=TRUE)) # dsCMatrix
str( dm9 <- dmperm(S9) )
(S9p <- with(dm9, S9[p, q]))
## looks good, but *not* quite upper triangular; these, too:
str( dm9.0 <- dmperm(S9, seed=-1)) # non-random too.
str( dm9_1 <- dmperm(S9, seed= 1)) # a random one
## The last two permutations differ, but have the same effect!
(S9p0 <- with(dm9.0, S9[p, q])) # .. hmm ..
stopifnot(all.equal(S9p0, S9p))# same as as default, but different from the random one


set.seed(11)
(M <- triu(rsparsematrix(9,11, 1/4)))
dM <- dmperm(M); with(dM, M[p, q])
(Mp <- M[sample.int(nrow(M)), sample.int(ncol(M))])
dMp <- dmperm(Mp); with(dMp, Mp[p, q])


set.seed(7)
(n7 <- rsparsematrix(5, 12, nnz = 10, rand.x = NULL))
str( dm.7 <- dmperm(n7) )
stopifnot(exprs = {
  lengths(dm.7[1:2]) == dim(n7)
  identical(dm.7,      dmperm(as(n7, "dMatrix")))
  identical(dm.7[1:4], dmperm(n7, nAns=4))
  identical(dm.7[1:2], dmperm(n7, nAns=2))
})
%% FIXME: Check permutations / blocks etc -- does it work ???
}
