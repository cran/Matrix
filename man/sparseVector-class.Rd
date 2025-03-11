\name{sparseVector-class}
\title{Sparse Vector Classes}
%
\docType{class}
\keyword{classes}
\keyword{manip}
%
\alias{sparseVector-class}
\alias{nsparseVector-class}
\alias{lsparseVector-class}
\alias{isparseVector-class}
\alias{dsparseVector-class}
\alias{zsparseVector-class}
% sparse
\alias{!,sparseVector-method}
\alias{Arith,sparseVector,ddenseMatrix-method}
\alias{Arith,sparseVector,dgeMatrix-method}
\alias{Arith,sparseVector,sparseVector-method}
\alias{Logic,sparseVector,dMatrix-method}
\alias{Logic,sparseVector,lMatrix-method}
\alias{Logic,sparseVector,nMatrix-method}
\alias{Logic,sparseVector,sparseVector-method}
\alias{Math,sparseVector-method}
\alias{Math2,sparseVector-method}
\alias{Ops,ANY,sparseVector-method}
\alias{Ops,matrix,sparseVector-method}
\alias{Ops,sparseVector,ANY-method}
\alias{Ops,sparseVector,Matrix-method}
\alias{Ops,sparseVector,matrix-method}
\alias{Ops,sparseVector,sparseVector-method}
\alias{Ops,sparseVector,vector-method}
\alias{Ops,vector,sparseVector-method}
\alias{Summary,sparseVector-method}
\alias{as.array,sparseVector-method}
\alias{as.complex,sparseVector-method}
\alias{as.integer,sparseVector-method}
\alias{as.logical,sparseVector-method}
\alias{as.matrix,sparseVector-method}
\alias{as.numeric,sparseVector-method}
\alias{as.vector,sparseVector-method}
\alias{coerce,ANY,sparseVector-method}
\alias{coerce,matrix,sparseVector-method}
\alias{coerce,sparseVector,CsparseMatrix-method}
\alias{coerce,sparseVector,Matrix-method}
\alias{coerce,sparseVector,RsparseMatrix-method}
\alias{coerce,sparseVector,TsparseMatrix-method}
\alias{coerce,sparseVector,denseMatrix-method}
\alias{coerce,sparseVector,dsparseVector-method}
\alias{coerce,sparseVector,generalMatrix-method}
\alias{coerce,sparseVector,isparseVector-method}
\alias{coerce,sparseVector,lsparseVector-method}
\alias{coerce,sparseVector,nsparseVector-method}
\alias{coerce,sparseVector,sparseMatrix-method}
\alias{coerce,sparseVector,unpackedMatrix-method}
\alias{coerce,sparseVector,zsparseVector-method}
\alias{coerce,vector,dsparseVector-method}
\alias{coerce,vector,isparseVector-method}
\alias{coerce,vector,lsparseVector-method}
\alias{coerce,vector,nsparseVector-method}
\alias{coerce,vector,sparseVector-method}
\alias{coerce,vector,zsparseVector-method}
\alias{diff,sparseVector-method}
\alias{dim<-,sparseVector-method}
\alias{head,sparseVector-method}
\alias{initialize,sparseVector-method}
\alias{length,sparseVector-method}
\alias{log,sparseVector-method}
\alias{mean,sparseVector-method}
\alias{rep,sparseVector-method}
\alias{show,sparseVector-method}
\alias{sort,sparseVector-method}
\alias{t,sparseVector-method}
\alias{tail,sparseVector-method}
\alias{toeplitz,sparseVector-method}
\alias{zapsmall,sparseVector-method}
% nsparse
\alias{!,nsparseVector-method}
\alias{which,nsparseVector-method}
% lsparse
\alias{!,lsparseVector-method}
\alias{Logic,lsparseVector,lsparseVector-method}
\alias{which,lsparseVector-method}
% isparse
% dsparse
\alias{-,dsparseVector,missing-method}
\alias{Arith,dsparseVector,dsparseVector-method}
% zsparse
%
\alias{c.sparseVector}
%
\description{Sparse Vector Classes: The virtual mother class
  \code{"sparseVector"} has the five actual daughter classes
   \code{"dsparseVector"}, \code{"isparseVector"},
   \code{"lsparseVector"}, \code{"nsparseVector"}, and
   \code{"zsparseVector"}, where we've mainly implemented methods for
   the \code{d*}, \code{l*} and \code{n*} ones.
}
\section{Slots}{
  \describe{
    \item{\code{length}:}{class \code{"numeric"} - the \code{\link{length}}
     of the sparse vector.  Note that \code{"numeric"} can be
     considerably larger than the maximal \code{"integer"},
     \code{\link{.Machine}$integer.max}, on purpose.}

    \item{\code{i}:}{class \code{"numeric"} - the (1-based) indices of
      the non-zero entries.  Must \emph{not} be \code{NA} and strictly
      sorted increasingly.

      Note that \code{"integer"} is \dQuote{part of} \code{"numeric"},
      and can (and often will) be used for non-huge sparseVectors.}

    \item{\code{x}:}{(for all but \code{"nsparseVector"}):
      the non-zero entries.  This is of class \code{"numeric"} for class
      \code{"dsparseVector"},  \code{"logical"} for class
      \code{"lsparseVector"}, etc.}
  }
}
\section{Methods}{
  \describe{
    \item{length}{\code{signature(x = "sparseVector")}: simply extracts
      the \code{length} slot.}
    \item{show}{\code{signature(object = "sparseVector")}: The
      \code{\link{show}} method for sparse vectors prints
      \emph{\dQuote{structural}} zeroes as \code{"."} using the
      non-exported \code{prSpVector} function which allows further
      customization such as replacing \code{"."} by \code{" "} (blank).

      Note that \code{\link{options}(max.print)} will influence how many
      entries of large sparse vectors are printed at all.}
    \item{as.vector}{\code{signature(x = "sparseVector", mode = "character")}
      coerces sparse vectors to \dQuote{regular}, i.e., atomic vectors.
      This is the same as \code{as(x, "vector")}.}
    \item{as}{..: see \code{coerce} below}
    \item{coerce}{\code{signature(from = "sparseVector", to = "sparseMatrix")}, and}
    \item{coerce}{\code{signature(from = "sparseMatrix", to = "sparseVector")},
      etc: coercions to and from sparse matrices (\code{\linkS4class{sparseMatrix}}) are
      provided and work analogously as in standard \R, i.e., a vector is
      coerced to a 1-column matrix.}
    \item{dim<-}{\code{signature(x = "sparseVector", value = "integer")}
      coerces a sparse vector to a sparse Matrix, i.e., an object
      inheriting from \code{\linkS4class{sparseMatrix}}, of the
      appropriate dimension.}
    \item{head}{\code{signature(x = "sparseVector")}: as with \R's
      (package \pkg{util}) \code{\link{head}}, \code{head(x,n)} (for
      \eqn{n >= 1}) is equivalent to \code{x[1:n]}, but here can be much
      more efficient, see the example.}
    \item{tail}{\code{signature(x = "sparseVector")}: analogous to
      \code{\link{head}}, see above.}
    \item{toeplitz}{\code{signature(x = "sparseVector")}: as
      \code{\link[stats]{toeplitz}(x)}, produce the \eqn{n \times n}
      Toeplitz matrix from \code{x}, where \code{n = length(x)}.}
    \item{rep}{\code{signature(x = "sparseVector")} repeat \code{x},
      with the same argument list \code{(x, times, length.out, each,
	...)} as the default method for rep().}
    \item{which}{\code{signature(x = "nsparseVector")} and}
    \item{which}{\code{signature(x = "lsparseVector")} return the
      indices of the non-zero entries (which is trivial for sparse vectors).}
    \item{Ops}{\code{signature(e1 = "sparseVector", e2 = "*")}: define
      arithmetic, compare and logic operations, (see
      \code{\link[=S4groupGeneric]{Ops}}).}
    \item{Summary}{\code{signature(x = "sparseVector")}: define
      all the \code{\link[=S4groupGeneric]{Summary}} methods.}
    \item{is.na, is.finite, is.infinite}{\code{(x = "sparseVector")}, and}
    \item{is.na, is.finite, is.infinite}{\code{(x = "nsparseVector")}:
      return \code{\link{logical}} or \code{"nsparseVector"} of the same
      length as \code{x}, indicating if/where \code{x} is
      \code{\link{NA}} (or \code{NaN}), finite or infinite, entirely
      analogously to the corresponding base \R functions.}
    \item{zapsmall}{\code{signature(x = "sparseVectors")}: typically used for
      numeric sparse vector: \code{\link{round}()} entries
      such that (relatively) very small entries become zero exactly.}
  }
  \code{c.sparseVector()} is an S3 method for all
  \code{"sparseVector"}s, but automatic dispatch only happens for the
  first argument, so it is useful also as regular \R function, see the
  examples.
}
\seealso{
  \code{\link{sparseVector}()} for friendly construction of sparse
  vectors (apart from \code{as(*, "sparseVector")}).
}
%\author{Martin}
\examples{
\dontshow{ % for R_DEFAULT_PACKAGES=NULL
library(utils, pos = "package:base", verbose = FALSE)
}
getClass("sparseVector")
getClass("dsparseVector")

sx <- c(0,0,3, 3.2, 0,0,0,-3:1,0,0,2,0,0,5,0,0)
(ss <- as(sx, "sparseVector"))

ix <- as.integer(round(sx))
(is <- as(ix, "sparseVector")) ## an "isparseVector" (!)
(ns <- sparseVector(i= c(7, 3, 2), length = 10)) # "nsparseVector"
## rep() works too:
(ri <- rep(is, length.out= 25))

## Using `dim<-`  as in base R :
r <- ss
dim(r) <- c(4,5) # becomes a sparse Matrix:
r
## or coercion (as as.matrix() in base R):
as(ss, "Matrix")
stopifnot(all(ss == print(as(ss, "CsparseMatrix"))))

## currently has "non-structural" FALSE -- printing as ":"
(lis <- is & FALSE)
(nn <- is[is == 0]) # all "structural" FALSE

## NA-case
sN <- sx; sN[4] <- NA
(svN <- as(sN, "sparseVector"))

v <- as(c(0,0,3, 3.2, rep(0,9),-3,0,-1, rep(0,20),5,0),
         "sparseVector")
v <- rep(rep(v, 50), 5000)
set.seed(1); v[sample(v@i, 1e6)] <- 0
str(v)
% Formal class 'dsparseVector' [package "Matrix"] with 3 slots
%   ..@ x     : num [1:250000] 3.2 -1 -3 3 5 3.2 -3 3 -1 5 ...
%   ..@ length: int 9500000
%   ..@ i     : int [1:250000] 4 16 52 155 189 194 204 231 244 265 ...
system.time(for(i in 1:4) hv <- head(v, 1e6))
##   user  system elapsed
##  0.033   0.000   0.032
system.time(for(i in 1:4) h2 <- v[1:1e6])
##   user  system elapsed
##  1.317   0.000   1.319

stopifnot(identical(hv, h2),
          identical(is | FALSE, is != 0),
          validObject(svN), validObject(lis), as.logical(is.na(svN[4])),
          identical(is^2 > 0, is & TRUE),
          all(!lis), !any(lis), length(nn@i) == 0, !any(nn), all(!nn),
          sum(lis) == 0, !prod(lis), range(lis) == c(0,0))

## create and use the t(.) method:
t(x20 <- sparseVector(c(9,3:1), i=c(1:2,4,7), length=20))
(T20 <- toeplitz(x20))
stopifnot(is(T20, "symmetricMatrix"), is(T20, "sparseMatrix"),
          identical(unname(as.matrix(T20)),
                    toeplitz(as.vector(x20))))

## c() method for "sparseVector" - also available as regular function
(c1 <- c(x20, 0,0,0, -10*x20))
(c2 <- c(ns, is, FALSE))
(c3 <- c(ns, !ns, TRUE, NA, FALSE))
(c4 <- c(ns, rev(ns)))
## here, c() would produce a list {not dispatching to c.sparseVector()}
(c5 <- c.sparseVector(0,0, x20))

## checking (consistency)
.v <- as.vector
.s <- function(v) as(v, "sparseVector")
stopifnot(exprs = {
    all.equal(c1, .s(c(.v(x20), 0,0,0, -10*.v(x20))),      tol = 0)
    all.equal(c2, .s(c(.v(ns), .v(is), FALSE)),            tol = 0)
    all.equal(c3, .s(c(.v(ns), !.v(ns), TRUE, NA, FALSE)), tol = 0)
    all.equal(c4, .s(c(.v(ns), rev(.v(ns)))),              tol = 0,
              check.class = FALSE)
    all.equal(c5, .s(c(0,0, .v(x20))),                     tol = 0)
})
}
