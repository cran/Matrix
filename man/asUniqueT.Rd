\name{asUniqueT}
\title{Standardize a Sparse Matrix in Triplet Format}
%
\keyword{array}
\keyword{logic}
\keyword{manip}
\keyword{utilities}
%
\alias{anyDuplicatedT}
\alias{isUniqueT}
\alias{asUniqueT}
\alias{aggregateT}
%
\alias{uniqTsparse}
%
\description{
  Detect or standardize a \code{\linkS4class{TsparseMatrix}} with
  unsorted or duplicated \eqn{(i,j)} pairs.
}
\usage{
anyDuplicatedT(x, \dots)
isUniqueT(x, byrow = FALSE, isT = is(x, "TsparseMatrix"))
asUniqueT(x, byrow = FALSE, isT = is(x, "TsparseMatrix"))
aggregateT(x)
}
\arguments{
  \item{x}{an \R{} object.  \code{anyDuplicatedT} and \code{aggregateT}
    require \code{x} inheriting from \code{\linkS4class{TsparseMatrix}}.
    \code{asUniqueT} requires \code{x} inheriting from
    \code{\linkS4class{Matrix}} and coerces \code{x}
    to \code{\linkS4class{TsparseMatrix}} if necessary.}
  \item{\dots}{optional arguments passed to the default method for
    generic function \code{\link{anyDuplicated}}.}
  \item{byrow}{a logical indicating if \code{x} should be sorted
    by row then by column.}
  \item{isT}{a logical indicating if \code{x} inherits from virtual
    class \code{\linkS4class{TsparseMatrix}}.}
}
\value{
  \code{anyDuplicatedT(x)} returns the index of the first duplicated
  \eqn{(i,j)} pair in \code{x} (0 if there are no duplicated pairs).

  \code{isUniqueT(x)} returns \code{TRUE} if \code{x} is a
  \code{\linkS4class{TsparseMatrix}} with sorted, nonduplicated
  \eqn{(i,j)} pairs and \code{FALSE} otherwise.

  \code{asUniqueT(x)} returns the unique
  \code{\linkS4class{TsparseMatrix}} representation of \code{x} with
  sorted, nonduplicated \eqn{(i,j)} pairs.  Values corresponding to
  identical \eqn{(i,j)} pairs are aggregated by addition, where in the
  logical case \dQuote{addition} refers to logical OR.
  
  \code{aggregateT(x)} aggregates without sorting.
}
\seealso{Virtual class \code{\linkS4class{TsparseMatrix}}.}
\examples{
\dontshow{ % for R_DEFAULT_PACKAGES=NULL
library(stats, pos = "package:base", verbose = FALSE)
library(utils, pos = "package:base", verbose = FALSE)
}
example("dgTMatrix-class", echo=FALSE)
## -> 'T2'  with (i,j,x) slots of length 5 each
T2u <- asUniqueT(T2)
stopifnot(## They "are" the same (and print the same):
          all.equal(T2, T2u, tol=0),
          ## but not internally:
          anyDuplicatedT(T2)  == 2,
          anyDuplicatedT(T2u) == 0,
          length(T2 @x) == 5,
          length(T2u@x) == 3)

isUniqueT(T2 ) # FALSE
isUniqueT(T2u) # TRUE

T3 <- T2u
T3[1, c(1,3)] <- 10; T3[2, c(1,5)] <- 20
T3u <- asUniqueT(T3)
str(T3u) # sorted in 'j', and within j, sorted in i
stopifnot(isUniqueT(T3u))

## Logical l.TMatrix and n.TMatrix :
(L2 <- T2 > 0)
validObject(L2u <- asUniqueT(L2))
(N2 <- as(L2, "nMatrix"))
validObject(N2u <- asUniqueT(N2))
stopifnot(N2u@i == L2u@i, L2u@i == T2u@i,  N2@i == L2@i, L2@i == T2@i,
          N2u@j == L2u@j, L2u@j == T2u@j,  N2@j == L2@j, L2@j == T2@j)
# now with a nasty NA  [partly failed in Matrix 1.1-5]:
L.0N <- L.1N <- L2
L.0N@x[1:2] <- c(FALSE, NA)
L.1N@x[1:2] <- c(TRUE, NA)
validObject(L.0N)
validObject(L.1N)
(m.0N <- as.matrix(L.0N))
(m.1N <- as.matrix(L.1N))
stopifnot(identical(10L, which(is.na(m.0N))), !anyNA(m.1N))
symnum(m.0N)
symnum(m.1N)
}
