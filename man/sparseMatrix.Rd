\name{sparseMatrix}
\alias{sparseMatrix}
\title{General Sparse Matrix Construction from Nonzero Entries}
\description{
  User friendly construction of a compressed, column-oriented, sparse
  matrix, inheriting from \code{\link{class}}
  \code{\linkS4class{CsparseMatrix}} (or
  \code{\linkS4class{TsparseMatrix}} if \code{giveCsparse} is false),
  from locations (and values) of its non-zero entries.

  This is the recommended user interface rather than direct
  \code{\link{new}("***Matrix", ....)} calls.
}
\usage{
sparseMatrix(i = ep, j = ep, p, x, dims, dimnames,
             symmetric = FALSE, triangular = FALSE, index1 = TRUE,
             repr = "C", giveCsparse = (repr == "C"),
             check = TRUE, use.last.ij = FALSE)
}
\arguments{
  \item{i,j}{integer vectors of the same length specifying the locations
    (row and column indices) of the non-zero (or non-\code{TRUE})
    entries of the matrix.  Note that for \emph{repeated} pairs
    \eqn{(i_k,j_k)}, when \code{x} is not missing, the corresponding
    \eqn{x_k} are \emph{added}, in consistency with the definition of
    the \code{"\linkS4class{TsparseMatrix}"} class, unless
    \code{use.last.ij} is true, in which case only the \emph{last} of
    the corresponding \eqn{(i_k, j_k, x_k)} triplet is used.}
  \item{p}{numeric (integer valued) vector of pointers, one for each
    column (or row), to the initial (zero-based) index of elements in the
    column (or row).  Exactly one of \code{i}, \code{j} or \code{p} must
    be missing.}
  \item{x}{
    optional values of the matrix entries.  If specified, must be of
    the same length as \code{i} / \code{j}, or of
    length one where it will be recycled to full length.  If missing,
    the resulting matrix will be a 0/1 patter\bold{n} matrix, i.e.,
    extending class \code{\linkS4class{nsparseMatrix}}.
  }
  \item{dims}{optional, non-negative, integer, dimensions vector of
    length 2.  Defaults to \code{c(max(i), max(j))}.}
  \item{dimnames}{optional list of \code{\link{dimnames}}; if not
    specified, none, i.e., \code{\link{NULL}} ones, are used.}
  \item{symmetric}{logical indicating if the resulting matrix should be
    symmetric.  In that case, only the lower or upper triangle needs to
    be specified via \eqn{(i/j/p)}.}
  \item{triangular}{logical indicating if the resulting matrix should be
    triangular.  In that case, the lower or upper triangle needs to
    be specified via \eqn{(i/j/p)}.}
  \item{index1}{logical scalar.  If \code{TRUE}, the default, the index
    vectors \code{i} and/or \code{j} are 1-based, as is the convention
    in \R.  That is, counting of rows and columns starts at 1.  If
    \code{FALSE} the index vectors are 0-based so counting of rows and
    columns starts at 0; this corresponds to the internal representation.}
  \item{repr}{\code{\link{character}} string, one of \code{"C"},
    \code{"T"}, or \code{"R"}, specifying the sparse \emph{repr}esentation to
    be used for the result, i.e., one from the super classes
    \code{\linkS4class{CsparseMatrix}}, \code{\linkS4class{TsparseMatrix}}, or
    \code{\linkS4class{RsparseMatrix}}.}
  \item{giveCsparse}{(\bold{deprecated}, replaced with \code{repr}):
    logical indicating if the result should be a
    \code{\linkS4class{CsparseMatrix}} or a
    \code{\linkS4class{TsparseMatrix}}, where the default was \code{TRUE},
    and now is determined from \code{repr}; very often Csparse matrices are
    more efficient subsequently, but not always.}
  \item{check}{logical indicating if a validity check is performed; do
    not set to \code{FALSE} unless you know what you're doing!}
  \item{use.last.ij}{logical indicating if in the case of repeated,
    i.e., duplicated pairs \eqn{(i_k, j_k)} only the last one should be
    used.  The default, \code{FALSE}, corresponds to the
    \code{"\linkS4class{TsparseMatrix}"} definition.}
}
\value{
  A sparse matrix, by default (from \code{repr = "C"}) in compressed,
  column-oriented form, as an \R object inheriting from both
  \code{\linkS4class{CsparseMatrix}} and \code{\linkS4class{generalMatrix}}.
}
\details{
  Exactly one of the arguments \code{i}, \code{j} and \code{p} must be
  missing.

  In typical usage, \code{p} is missing, \code{i} and \code{j} are
  vectors of positive integers and \code{x} is a numeric vector.  These
  three vectors, which must have the same length, form the triplet
  representation of the sparse matrix.

  If \code{i} or \code{j} is missing then \code{p} must be a
  non-decreasing integer vector whose first element is zero.  It
  provides the compressed, or \dQuote{pointer} representation of the row
  or column indices, whichever is missing.  The expanded form of \code{p},
  \code{rep(seq_along(dp),dp)} where \code{dp <- diff(p)}, is used as
  the (1-based) row or column indices.

  You cannot set both \code{singular} and \code{triangular} to true;
  rather use \code{\link{Diagonal}()} (or its alternatives, see there).

  The values of \code{i}, \code{j}, \code{p} and \code{index1} are used
  to create 1-based index vectors \code{i} and \code{j} from which a
  \code{\linkS4class{TsparseMatrix}} is constructed, with numerical
  values given by \code{x}, if non-missing.  Note that in that case,
  when some pairs \eqn{(i_k,j_k)} are repeated (aka
  \dQuote{duplicated}), the corresponding \eqn{x_k} are \emph{added}, in
  consistency with the definition of the
  \code{"\linkS4class{TsparseMatrix}"} class, unless \code{use.last.ij}
  is set to true.
  %%
  By default, when \code{repr = "C"}, the \code{\linkS4class{CsparseMatrix}}
  derived from this triplet form is returned, where \code{repr = "R"} now
  allows to directly get an \code{\linkS4class{RsparseMatrix}} and
  \code{repr = "T"} leaves the result as \code{\linkS4class{TsparseMatrix}}.

  The reason for returning a \code{\linkS4class{CsparseMatrix}} object
  instead of the triplet format by default is that the compressed column
  form is easier to work with when performing matrix operations.  In
  particular, if there are no zeros in \code{x} then a
  \code{\linkS4class{CsparseMatrix}} is a unique representation of the
  sparse matrix.
}
\note{% We say so above (-> {index1}), but some do not read that
  You \emph{do} need to use \code{index1 = FALSE} (or add \code{+ 1}
  to \code{i} and \code{j}) if you want use the 0-based \code{i} (and
  \code{j}) slots from existing sparse matrices.
}
\seealso{\code{\link{Matrix}(*, sparse=TRUE)} for the constructor of
  such matrices from a \emph{dense} matrix.  That is easier in small
  sample, but much less efficient (or impossible) for large matrices,
  where something like \code{sparseMatrix()} is needed.
  Further \code{\link{bdiag}} and \code{\link{Diagonal}} for (block-)diagonal and
  \code{\link{bandSparse}} for banded sparse matrix constructors.

  Random sparse matrices via \code{\link{rsparsematrix}()}.

  The standard \R \code{\link{xtabs}(*, sparse=TRUE)}, for sparse tables
  and \code{\link{sparse.model.matrix}()} for building sparse model
  matrices.

  Consider \code{\linkS4class{CsparseMatrix}} and similar class
  definition help files.
}
\examples{
## simple example
i <- c(1,3:8); j <- c(2,9,6:10); x <- 7 * (1:7)
(A <- sparseMatrix(i, j, x = x))                    ##  8 x 10 "dgCMatrix"
summary(A)
str(A) # note that *internally* 0-based row indices are used

(sA <- sparseMatrix(i, j, x = x, symmetric = TRUE)) ## 10 x 10 "dsCMatrix"
(tA <- sparseMatrix(i, j, x = x, triangular= TRUE)) ## 10 x 10 "dtCMatrix"
stopifnot( all(sA == tA + t(tA)) ,
           identical(sA, as(tA + t(tA), "symmetricMatrix")))

## dims can be larger than the maximum row or column indices
(AA <- sparseMatrix(c(1,3:8), c(2,9,6:10), x = 7 * (1:7), dims = c(10,20)))
summary(AA)

## i, j and x can be in an arbitrary order, as long as they are consistent
set.seed(1); (perm <- sample(1:7))
(A1 <- sparseMatrix(i[perm], j[perm], x = x[perm]))
stopifnot(identical(A, A1))

## The slots are 0-index based, so
try( sparseMatrix(i=A@i, p=A@p, x= seq_along(A@x)) )
## fails and you should say so: 1-indexing is FALSE:
     sparseMatrix(i=A@i, p=A@p, x= seq_along(A@x), index1 = FALSE)

## the (i,j) pairs can be repeated, in which case the x's are summed
(args <- data.frame(i = c(i, 1), j = c(j, 2), x = c(x, 2)))
(Aa <- do.call(sparseMatrix, args))
## explicitly ask for elimination of such duplicates, so
## that the last one is used:
(A. <- do.call(sparseMatrix, c(args, list(use.last.ij = TRUE))))
stopifnot(Aa[1,2] == 9, # 2+7 == 9
          A.[1,2] == 2) # 2 was *after* 7

## for a pattern matrix, of course there is no "summing":
(nA <- do.call(sparseMatrix, args[c("i","j")]))

dn <- list(LETTERS[1:3], letters[1:5])
## pointer vectors can be used, and the (i,x) slots are sorted if necessary:
m <- sparseMatrix(i = c(3,1, 3:2, 2:1), p= c(0:2, 4,4,6), x = 1:6, dimnames = dn)
m
str(m)
stopifnot(identical(dimnames(m), dn))

sparseMatrix(x = 2.72, i=1:3, j=2:4) # recycling x
sparseMatrix(x = TRUE, i=1:3, j=2:4) # recycling x, |--> "lgCMatrix"

## no 'x' --> patter*n* matrix:
(n <- sparseMatrix(i=1:6, j=rev(2:7)))# -> ngCMatrix

## an empty sparse matrix:
(e <- sparseMatrix(dims = c(4,6), i={}, j={}))

## a symmetric one:
(sy <- sparseMatrix(i= c(2,4,3:5), j= c(4,7:5,5), x = 1:5,
                    dims = c(7,7), symmetric=TRUE))
stopifnot(isSymmetric(sy),
          identical(sy, ## switch i <-> j {and transpose }
    t( sparseMatrix(j= c(2,4,3:5), i= c(4,7:5,5), x = 1:5,
                    dims = c(7,7), symmetric=TRUE))))

## rsparsematrix() calls sparseMatrix() :
M1 <- rsparsematrix(1000, 20, nnz = 200)
summary(M1)

## pointers example in converting from other sparse matrix representations.
if(require(SparseM) && packageVersion("SparseM") >= 0.87 &&
   nzchar(dfil <- system.file("extdata", "rua_32_ax.rua", package = "SparseM"))) {
  X <- model.matrix(read.matrix.hb(dfil))
  XX <- sparseMatrix(j = X@ja, p = X@ia - 1L, x = X@ra, dims = X@dimension)
  validObject(XX)

  ## Alternatively, and even more user friendly :
  X. <- as(X, "Matrix")  # or also
  X2 <- as(X, "sparseMatrix")
  stopifnot(identical(XX, X.), identical(X., X2))
}% if
}% example
\keyword{array}
