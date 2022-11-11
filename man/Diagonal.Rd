\name{Diagonal}
\title{Construct a Diagonal Matrix}
%
\alias{Diagonal}
\alias{.sparseDiagonal}
\alias{.trDiagonal}
\alias{.symDiagonal}
%
\description{
  Construct a formally diagonal \code{\linkS4class{Matrix}},
  i.e., an object inheriting from virtual class
  \code{\linkS4class{diagonalMatrix}}
  (or, if desired, a \emph{mathematically} diagonal
  \code{\linkS4class{CsparseMatrix}}).
}
\usage{
Diagonal(n, x = NULL, names = FALSE)

.sparseDiagonal(n, x = NULL, uplo = "U", shape = "t", unitri = TRUE, kind, cols)
    .trDiagonal(n, x = NULL, uplo = "U", unitri = TRUE, kind)
   .symDiagonal(n, x = NULL, uplo = "U", kind)
}
\arguments{
  \item{n}{integer indicating the dimension of the (square) matrix.
    If missing, then \code{length(x)} is used.}
  \item{x}{numeric or logical vector listing values for the diagonal
    entries, to be recycled as necessary.  If \code{NULL} (the default),
    then the result is a unit diagonal matrix.  \code{.sparseDiagonal()}
    and friends ignore non-\code{NULL} \code{x} when \code{kind = "n"}.}
  \item{names}{either \code{\link{logical}} \code{TRUE} or \code{FALSE} or
    then a \code{\link{character}} vector of \code{\link{length}}
    \code{n}.  If true \emph{and} \code{\link{names}(x)} is not
    \code{NULL}, use that as both row and column names for the resulting
    matrix.  When a character vector, use it for both dimnames.}
  \item{uplo}{one of \code{c("U","L")}, specifying the \code{uplo} slot
    of the result if the result is formally triangular of symmetric.}
  \item{shape}{one of \code{c("t","s","g")}, indicating if the result
    should be formally triangular, symmetric, or \dQuote{general}.
    The result will inherit from virtual class
    \code{\linkS4class{triangularMatrix}},
    \code{\linkS4class{symmetricMatrix}}, or
    \code{\linkS4class{generalMatrix}}, respectively.}
  \item{unitri}{logical indicating if a formally triangular result with
    ones on the diagonal should be formally \emph{unit} triangular, i.e.,
    with \code{diag} slot equal to \code{"U"} rather than \code{"N"}.}
  \item{kind}{one of \code{c("d","l","n")}, indicating the \dQuote{mode}
    of the result: numeric, logical, or pattern.
    The result will inherit from virtual class
    \code{\linkS4class{dsparseMatrix}},
    \code{\linkS4class{lsparseMatrix}}, or
    \code{\linkS4class{nsparseMatrix}}, respectively.
    Values other than \code{"n"} are ignored when \code{x} is
    non-\code{NULL}; in that case the mode is determined by
    \code{\link{typeof}(x)}.}
  \item{cols}{optional integer vector with values in \code{0:(n-1)},
    indexing columns of the specified diagonal matrix.  If specified,
    then the result is (mathematically) \code{D[, cols+1]} rather
    than \code{D}, where \code{D = Diagonal(n, x)}, and it is always
    \dQuote{general} (i.e., \code{shape} is ignored).}
}
\value{
  \code{Diagonal()} returns an object inheriting from virtual class
  \code{\linkS4class{diagonalMatrix}}.

  \code{.sparseDiagonal()} returns a \code{\linkS4class{CsparseMatrix}}
  representation of \code{Diagonal(n, x)} or, if \code{cols} is given,
  of \code{Diagonal(n, x)[, cols+1]}.  The precise class of the result
  depends on \code{shape} and \code{kind}.

  \code{.trDiagonal()} and \code{.symDiagonal()} are simple wrappers,
  for \code{.sparseDiagonal(shape = "t")} and
  \code{.sparseDiagonal(shape = "s")}, respectively.

  \code{.sparseDiagonal()} exists primarily to leverage efficient
  C-level methods available for \code{CsparseMatrix}.
}
\author{Martin Maechler}
\seealso{the generic function \code{\link{diag}} for \emph{extraction}
  of the diagonal from a matrix works for all \dQuote{Matrices}.

  \code{\link{bandSparse}} constructs a \emph{banded} sparse matrix from
  its non-zero sub-/super - diagonals.  \code{\link{band}(A)} returns a
  band matrix containing some sub-/super - diagonals of \code{A}.

  \code{\link{Matrix}} for general matrix construction;
  further, class \code{\linkS4class{diagonalMatrix}}.
}
\examples{
Diagonal(3)
Diagonal(x = 10^(3:1))
Diagonal(x = (1:4) >= 2)#-> "ldiMatrix"

## Use Diagonal() + kronecker() for "repeated-block" matrices:
M1 <- Matrix(0+0:5, 2,3)
(M <- kronecker(Diagonal(3), M1))

(S <- crossprod(Matrix(rbinom(60, size=1, prob=0.1), 10,6)))
(SI <- S + 10*.symDiagonal(6)) # sparse symmetric still
stopifnot(is(SI, "dsCMatrix"))
(I4 <- .sparseDiagonal(4, shape="t"))# now (2012-10) unitriangular
stopifnot(I4@diag == "U", all(I4 == diag(4)))
\dontshow{% checking some "unit-diagonality":
  L <- Diagonal(5, TRUE)
  stopifnot(L@diag == "U", identical(L, Diagonal(5) > 0))
}
}
\keyword{array}
\keyword{algebra}
