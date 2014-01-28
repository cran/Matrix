\name{dgTMatrix-class}
\title{Sparse matrices in triplet form}
\docType{class}
\alias{dgTMatrix-class}
\alias{+,dgTMatrix,dgTMatrix-method}
\alias{coerce,dgTMatrix,dgCMatrix-method}
\alias{coerce,dgTMatrix,dgeMatrix-method}
\alias{coerce,dgTMatrix,matrix-method}
\alias{coerce,dgTMatrix,dtCMatrix-method}
\alias{coerce,dgTMatrix,dsTMatrix-method}
\alias{coerce,dgTMatrix,dtTMatrix-method}
\alias{coerce,dgTMatrix,symmetricMatrix-method}
\alias{coerce,dgTMatrix,triangularMatrix-method}
\alias{coerce,dgeMatrix,dgTMatrix-method}
\alias{coerce,matrix,dgTMatrix-method}
\description{The \code{"dgTMatrix"} class is the class of sparse
  matrices stored as (possibly redundant) triplets.  The internal
  representation is not at all unique, contrary to the one for class
  \code{\linkS4class{dgCMatrix}}.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("dgTMatrix",
    ...)}, but more typically via \code{as(*, "dgTMatrix")},
  \code{\link{spMatrix}()}, or
  \code{\link{sparseMatrix}(*, giveCsparse=FALSE)}.
}
\section{Slots}{
  \describe{
    \item{\code{i}:}{Object of class \code{"integer"} - the row indices
      of non-zero entries \emph{in 0-base}, i.e., must be in
      \code{0:(nrow(.)-1)}.}
    \item{\code{j}:}{Object of class \code{"integer"} - the column
      indices of non-zero entries.  Must be the same length as slot
      \code{i} and \emph{0-based} as well, i.e., in \code{0:(ncol(.)-1)}.}
    \item{\code{x}:}{Object of class \code{"numeric"} - the (non-zero)
      entry at position \code{(i,j)}.  Must be the same length as slot
      \code{i}. If an index pair occurs more than once the corresponding
      values of slot \code{x} are added to form the element of the matrix.}
    \item{\code{Dim}:}{Object of class \code{"integer"} of length 2 -
      the dimensions of the matrix.}
  }
}
\section{Methods}{
  \describe{
    \item{+}{\code{signature(e1 = "dgTMatrix", e2 = "dgTMatrix")}}
    \item{coerce}{\code{signature(from = "dgTMatrix", to = "dgCMatrix")}}
    \item{coerce}{\code{signature(from = "dgTMatrix", to = "dgeMatrix")}}
    \item{coerce}{\code{signature(from = "dgTMatrix", to = "matrix")},
      and typically coercion methods for more specific signatures, we
      are not mentioning here.

      Note that these are not guaranteed to continue to exist,
      but rather you should use calls like \code{as(x,
	"CsparseMatrix")}, \code{as(x, "generalMatrix")},
      \code{as(x, "dMatrix")}, i.e. coercion to higher level virtual classes.}
    \item{coerce}{\code{signature(from = "matrix", to = "dgTMatrix")},
      (direct coercion from tradition matrix).}

    \item{image}{\code{signature(x = "dgTMatrix")}: plots an image of
      \code{x} using the \code{\link[lattice]{levelplot}} function}
    \item{t}{\code{signature(x = "dgTMatrix")}: returns the transpose of
      \code{x}}
  }
}
%\references{}
%\author{}
\note{Triplet matrices are a convenient form in which to construct sparse
  matrices after which they can be coerced to
  \code{\linkS4class{dgCMatrix}} objects.

  Note that both \code{new(.)} and \code{\link{spMatrix}} constructors
  for \code{"dgTMatrix"} (and other \code{"linkS4class{TsparseMatrix}"}
  classes) implicitly add \eqn{x_k}'s that belong to identical
  \eqn{(i_k, j_k)} pairs.
}
\seealso{
  Class \code{\linkS4class{dgCMatrix}} or the superclasses
  \code{\linkS4class{dsparseMatrix}} and
  \code{\linkS4class{TsparseMatrix}};
}
\examples{
m <- Matrix(0+1:28, nrow = 4)
m[-3,c(2,4:5,7)] <- m[ 3, 1:4] <- m[1:3, 6] <- 0
(mT <- as(m, "dgTMatrix"))
str(mT)
mT[1,]
mT[4, drop = FALSE]
stopifnot(identical(mT[lower.tri(mT)],
                    m [lower.tri(m) ]))
mT[lower.tri(mT,diag=TRUE)] <- 0
mT

## Triplet representation with repeated (i,j) entries
## *adds* the corresponding x's:
T2 <- new("dgTMatrix",
          i = as.integer(c(1,1,0,3,3)),
          j = as.integer(c(2,2,4,0,0)), x=10*1:5, Dim=4:5)
str(T2) # contains (i,j,x) slots exactly as above, but
T2 ## has only three non-zero entries, as for repeated (i,j)'s,
   ## the corresponding x's are "implicitly" added
stopifnot(nnzero(T2) == 3)
}
\keyword{classes}
\keyword{algebra}
