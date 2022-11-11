\name{sparseMatrix}
\title{General Sparse Matrix Construction from Nonzero Entries}
\alias{sparseMatrix}
\description{
  User-friendly construction of sparse matrices (inheriting from
  virtual \code{\link{class}} \code{\linkS4class{CsparseMatrix}},
  \code{\linkS4class{RsparseMatrix}}, or
  \code{\linkS4class{TsparseMatrix}})
  from the positions and values of their nonzero entries.

  This interface is recommended over direct construction via
  calls such as \code{\link{new}("..[CRT]Matrix", ...)}.
}
\usage{
sparseMatrix(i, j, p, x, dims, dimnames,
             symmetric = FALSE, triangular = FALSE, index1 = TRUE,
             repr = c("C", "R", "T"), giveCsparse,
             check = TRUE, use.last.ij = FALSE)
}
\arguments{
  \item{i,j}{integer vectors of equal length specifying the positions
    (row and column indices) of the nonzero (or non-\code{TRUE}) entries
    of the matrix.  Note that, when \code{x} is non-missing, the
    \eqn{x_k} corresponding to \emph{repeated} pairs \eqn{(i_k,j_k)}
    are \emph{added}, for consistency with the definition of class
    \code{\linkS4class{TsparseMatrix}}, unless \code{use.last.ij} is
    \code{TRUE}, in which case only the \emph{last} such \eqn{x_k} is
    used.}
  \item{p}{integer vector of pointers, one for each column (or row),
    to the initial (zero-based) index of elements in the column (or row).
    Exactly one of \code{i}, \code{j}, and \code{p} must be missing.}
  \item{x}{optional, typically nonzero values for the matrix entries.
    If specified, then the length must equal that of \code{i}
    (or \code{j}) or equal 1, in which case \code{x} is recycled as
    necessary.  If missing, then the result is a \bold{n}onzero pattern
    matrix, i.e., inheriting from class \code{\linkS4class{nsparseMatrix}}.}
  \item{dims}{optional length-2 integer vector of matrix dimensions.
    If missing, then \code{!index1+c(max(i),max(j))} is used.}
  \item{dimnames}{optional list of \code{\link{dimnames}}; if missing,
    then \code{\link{NULL}} ones are used.}
  \item{symmetric}{logical indicating if the resulting matrix should
    be symmetric.  In that case, \eqn{(i,j,p)} should specify only one
    triangle (upper or lower).}
  \item{triangular}{logical indicating if the resulting matrix should
    be triangular.  In that case, \eqn{(i,j,p)} should specify only one
    triangle (upper or lower).}
  \item{index1}{logical.  If \code{TRUE} (the default), then \code{i}
    and \code{j} are interpreted as 1-based indices, following the \R
    convention.  That is, counting of rows and columns starts at 1.
    If \code{FALSE}, then they are interpreted as 0-based indices.}
  \item{repr}{\code{\link{character}} string, one of \code{"C"},
    \code{"R"}, and \code{"T"}, specifying the \bold{repr}esentation
    of the sparse matrix result, i.e., specifying one of the virtual
    classes \code{\linkS4class{CsparseMatrix}},
    \code{\linkS4class{RsparseMatrix}}, and
    \code{\linkS4class{TsparseMatrix}}.}
  \item{giveCsparse}{(\bold{deprecated}, replaced by \code{repr})
    logical indicating if the result should inherit from
    \code{\linkS4class{CsparseMatrix}} or
    \code{\linkS4class{TsparseMatrix}}.
    Note that operations involving \code{CsparseMatrix} are very often
    (but not always) more efficient.}
  \item{check}{logical indicating whether to check that the result is
    formally valid before returning.  Do not set to \code{FALSE} unless
    you know what you are doing!}
  \item{use.last.ij}{logical indicating if, in the case of repeated
    (duplicated) pairs \eqn{(i_k,j_k)}, only the last pair should be
    used.  \code{FALSE} (the default) is consistent with the definiton
    of class \code{\linkS4class{TsparseMatrix}}.}
}
\value{
  A sparse matrix, by default in compressed sparse column format and
  (formally) without symmetric or triangular structure, i.e.,
  by default inheriting from both \code{\linkS4class{CsparseMatrix}}
  and \code{\linkS4class{generalMatrix}}.
}
\details{ % FIXME: quite a bit of repetition here (and below) ...
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
  \code{\linkS4class{TsparseMatrix}} class, unless \code{use.last.ij}
  is set to true.
  
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
