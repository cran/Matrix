\name{nsparseMatrix-classes}
\title{Sparse "pattern" Matrices}
\docType{class}
%
\alias{nsparseMatrix-class}
\alias{ngCMatrix-class}
\alias{ngRMatrix-class}
\alias{ngTMatrix-class}
\alias{ntCMatrix-class}
\alias{ntRMatrix-class}
\alias{ntTMatrix-class}
\alias{nsCMatrix-class}
\alias{nsRMatrix-class}
\alias{nsTMatrix-class}
% nsparse
\alias{!,nsparseMatrix-method}
\alias{-,nsparseMatrix,missing-method}
\alias{Arith,nsparseMatrix,Matrix-method}
\alias{Arith,dsparseMatrix,nsparseMatrix-method}
\alias{Arith,lsparseMatrix,nsparseMatrix-method}
\alias{Arith,nsparseMatrix,dsparseMatrix-method}
\alias{Arith,nsparseMatrix,lsparseMatrix-method}
\alias{Ops,nsparseMatrix,dsparseMatrix-method}
\alias{Ops,nsparseMatrix,lsparseMatrix-method}
\alias{Ops,nsparseMatrix,sparseMatrix-method}
\alias{as.logical,nsparseMatrix-method}
\alias{as.numeric,nsparseMatrix-method}
\alias{coerce,matrix,nsparseMatrix-method}
\alias{coerce,nsparseMatrix,dMatrix-method}
\alias{coerce,nsparseMatrix,dsparseMatrix-method}
\alias{coerce,nsparseMatrix,indMatrix-method}
\alias{coerce,nsparseMatrix,lMatrix-method}
\alias{coerce,nsparseMatrix,lsparseMatrix-method}
\alias{coerce,nsparseMatrix,pMatrix-method}
\alias{coerce,numLike,nsparseMatrix-method}
\alias{which,nsparseMatrix-method}
% ngC
\alias{coerce,matrix,ngCMatrix-method}
\alias{coerce,ngCMatrix,dgCMatrix-method}
\alias{coerce,ngCMatrix,lgCMatrix-method}
\alias{coerce,ngCMatrix,ngRMatrix-method}
\alias{coerce,ngCMatrix,ngTMatrix-method}
\alias{coerce,ngCMatrix,ngeMatrix-method}
\alias{coerce,ngCMatrix,nsCMatrix-method}
\alias{coerce,ngCMatrix,nsTMatrix-method}
\alias{coerce,ngCMatrix,ntCMatrix-method}
\alias{coerce,ngCMatrix,ntTMatrix-method}
\alias{coerce,numLike,ngCMatrix-method}
% ngR
\alias{coerce,matrix,ngRMatrix-method}
\alias{coerce,ngRMatrix,dgRMatrix-method}
\alias{coerce,ngRMatrix,lgRMatrix-method}
\alias{coerce,ngRMatrix,ngCMatrix-method}
\alias{coerce,ngRMatrix,ngTMatrix-method}
\alias{coerce,ngRMatrix,ngeMatrix-method}
\alias{coerce,ngRMatrix,nsRMatrix-method}
\alias{coerce,ngRMatrix,ntRMatrix-method}
\alias{coerce,numLike,ngRMatrix-method}
% ngT
\alias{coerce,matrix,ngTMatrix-method}
\alias{coerce,ngTMatrix,dgTMatrix-method}
\alias{coerce,ngTMatrix,lgTMatrix-method}
\alias{coerce,ngTMatrix,lgeMatrix-method}
\alias{coerce,ngTMatrix,ngCMatrix-method}
\alias{coerce,ngTMatrix,ngRMatrix-method}
\alias{coerce,ngTMatrix,ngeMatrix-method}
\alias{coerce,ngTMatrix,nsCMatrix-method}
\alias{coerce,ngTMatrix,nsTMatrix-method}
\alias{coerce,ngTMatrix,ntCMatrix-method}
\alias{coerce,ngTMatrix,ntTMatrix-method}
\alias{coerce,numLike,ngTMatrix-method}
\alias{which,ngTMatrix-method}
% ntC
\alias{coerce,matrix,ntCMatrix-method}
\alias{coerce,ntCMatrix,dtCMatrix-method}
\alias{coerce,ntCMatrix,ltCMatrix-method}
\alias{coerce,ntCMatrix,ngCMatrix-method}
\alias{coerce,ntCMatrix,ngTMatrix-method}
\alias{coerce,ntCMatrix,nsCMatrix-method}
\alias{coerce,ntCMatrix,ntRMatrix-method}
\alias{coerce,ntCMatrix,ntTMatrix-method}
\alias{coerce,ntCMatrix,ntpMatrix-method}
\alias{coerce,ntCMatrix,ntrMatrix-method}
% ntR
\alias{coerce,matrix,ntRMatrix-method}
\alias{coerce,ntRMatrix,dtRMatrix-method}
\alias{coerce,ntRMatrix,ltRMatrix-method}
\alias{coerce,ntRMatrix,ngRMatrix-method}
\alias{coerce,ntRMatrix,nsRMatrix-method}
\alias{coerce,ntRMatrix,ntCMatrix-method}
\alias{coerce,ntRMatrix,ntTMatrix-method}
\alias{coerce,ntRMatrix,ntpMatrix-method}
\alias{coerce,ntRMatrix,ntrMatrix-method}
% ntT
\alias{coerce,matrix,ntTMatrix-method}
\alias{coerce,ntTMatrix,dtTMatrix-method}
\alias{coerce,ntTMatrix,ltTMatrix-method}
\alias{coerce,ntTMatrix,ngCMatrix-method}
\alias{coerce,ntTMatrix,ngTMatrix-method}
\alias{coerce,ntTMatrix,nsTMatrix-method}
\alias{coerce,ntTMatrix,ntCMatrix-method}
\alias{coerce,ntTMatrix,ntRMatrix-method}
\alias{coerce,ntTMatrix,ntpMatrix-method}
\alias{coerce,ntTMatrix,ntrMatrix-method}
\alias{which,ntTMatrix-method}
% nsC
\alias{coerce,matrix,nsCMatrix-method}
\alias{coerce,nsCMatrix,RsparseMatrix-method}
\alias{coerce,nsCMatrix,dsCMatrix-method}
\alias{coerce,nsCMatrix,lsCMatrix-method}
\alias{coerce,nsCMatrix,ngCMatrix-method}
\alias{coerce,nsCMatrix,ngTMatrix-method}
\alias{coerce,nsCMatrix,nsRMatrix-method}
\alias{coerce,nsCMatrix,nsTMatrix-method}
\alias{coerce,nsCMatrix,nspMatrix-method}
\alias{coerce,nsCMatrix,nsyMatrix-method}
\alias{coerce,nsCMatrix,ntCMatrix-method}
% nsR
\alias{coerce,matrix,nsRMatrix-method}
\alias{coerce,nsRMatrix,CsparseMatrix-method}
\alias{coerce,nsRMatrix,dsRMatrix-method}
\alias{coerce,nsRMatrix,lsRMatrix-method}
\alias{coerce,nsRMatrix,ngRMatrix-method}
\alias{coerce,nsRMatrix,nsCMatrix-method}
\alias{coerce,nsRMatrix,nsTMatrix-method}
\alias{coerce,nsRMatrix,nspMatrix-method}
\alias{coerce,nsRMatrix,nsyMatrix-method}
\alias{coerce,nsRMatrix,ntRMatrix-method}
% nsT
\alias{coerce,matrix,nsTMatrix-method}
\alias{coerce,nsTMatrix,dsTMatrix-method}
\alias{coerce,nsTMatrix,lsTMatrix-method}
\alias{coerce,nsTMatrix,ngCMatrix-method}
\alias{coerce,nsTMatrix,ngTMatrix-method}
\alias{coerce,nsTMatrix,nsCMatrix-method}
\alias{coerce,nsTMatrix,nsRMatrix-method}
\alias{coerce,nsTMatrix,nspMatrix-method}
\alias{coerce,nsTMatrix,nsyMatrix-method}
\alias{coerce,nsTMatrix,ntTMatrix-method}
\alias{which,nsTMatrix-method}
%
\description{The \code{nsparseMatrix} class is a virtual class of sparse
  \emph{\dQuote{pattern}} matrices, i.e., binary matrices conceptually
  with \code{TRUE}/\code{FALSE} entries.  Only the positions of the
  elements that are \code{TRUE} are stored.

  These can be stored in the \dQuote{triplet} form
  (\code{\linkS4class{TsparseMatrix}}, subclasses \code{ngTMatrix},
  \code{nsTMatrix}, and \code{ntTMatrix} which really contain pairs, not
  triplets) or in compressed column-oriented form (class
  \code{\linkS4class{CsparseMatrix}}, subclasses \code{ngCMatrix},
  \code{nsCMatrix}, and \code{ntCMatrix}) or--\emph{rarely}--in
  compressed row-oriented form (class \code{\linkS4class{RsparseMatrix}},
  subclasses \code{ngRMatrix}, \code{nsRMatrix}, and \code{ntRMatrix}).
  The second letter in the name of these non-virtual classes indicates
  \code{g}eneral, \code{s}ymmetric, or \code{t}riangular.

}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("ngCMatrix",
    ...)} and so on.  More frequently objects are created by coercion of
  a numeric sparse matrix to the pattern form for use in
  the symbolic analysis phase
  of an algorithm involving sparse matrices.  Such algorithms often
  involve two phases: a symbolic phase wherein the positions of the
  non-zeros in the result are determined and a numeric phase wherein the
  actual results are calculated.  During the symbolic phase only the
  positions of the non-zero elements in any operands are of interest,
  hence numeric sparse matrices can be treated as sparse pattern
  matrices.
}
\section{Slots}{
  \describe{
    \item{\code{uplo}:}{Object of class \code{"character"}. Must be
      either "U", for upper triangular, and "L", for lower
      triangular. Present in the triangular and symmetric classes but not
      in the general class.}
    \item{\code{diag}:}{Object of class \code{"character"}. Must be
      either \code{"U"}, for unit triangular (diagonal is all ones), or
      \code{"N"} for non-unit.  The implicit diagonal elements are not
      explicitly stored when \code{diag} is \code{"U"}.  Present in the
      triangular classes only.}
    \item{\code{p}:}{Object of class \code{"integer"} of pointers, one
      for each column (row), to the initial (zero-based) index of elements in
      the column.  Present in compressed column-oriented and compressed
      row-oriented forms only.}
    \item{\code{i}:}{Object of class \code{"integer"} of length nnzero
      (number of non-zero elements).  These are the row numbers for
      each TRUE element in the matrix.  All other elements are FALSE.
      Present in triplet and compressed column-oriented forms only.}
    \item{\code{j}:}{Object of class \code{"integer"} of length nnzero
      (number of non-zero elements).  These are the column numbers for
      each TRUE element in the matrix.  All other elements are FALSE.
      Present in triplet and compressed row-oriented forms only.}
    \item{\code{Dim}:}{Object of class \code{"integer"} - the dimensions
      of the matrix.}
  }
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "dgCMatrix", to =
	"ngCMatrix")}, and many similar ones; typically you should
      coerce to \code{"nsparseMatrix"} (or \code{"nMatrix"}).  Note that
      coercion to a sparse pattern matrix records all the potential
      non-zero entries, i.e., explicit (\dQuote{non-structural}) zeroes
      are coerced to \code{TRUE}, not \code{FALSE}, see the example.
    }
    \item{t}{\code{signature(x = "ngCMatrix")}: returns the transpose
      of \code{x}}

    \item{which}{\code{signature(x = "lsparseMatrix")}, semantically
      equivalent to \pkg{base} function \code{\link{which}(x, arr.ind)};
      for details, see the \code{\linkS4class{lMatrix}} class documentation.}
  }
}
%\references{}
%\author{}
%\note{}
\seealso{
  the class \code{\linkS4class{dgCMatrix}}
}
\examples{
(m <- Matrix(c(0,0,2:0), 3,5, dimnames=list(LETTERS[1:3],NULL)))
## ``extract the nonzero-pattern of (m) into an nMatrix'':
nm <- as(m, "nsparseMatrix") ## -> will be a "ngCMatrix"
str(nm) # no 'x' slot
nnm <- !nm # no longer sparse
## consistency check:
stopifnot(xor(as( nm, "matrix"),
              as(nnm, "matrix")))

## low-level way of adding "non-structural zeros" :
nnm <- as(nnm, "lsparseMatrix") # "lgCMatrix"
nnm@x[2:4] <- c(FALSE, NA, NA)
nnm
as(nnm, "nMatrix") # NAs *and* non-structural 0  |--->  'TRUE'

data(KNex)
nmm <- as(KNex $ mm, "nMatrix")
str(xlx <- crossprod(nmm))# "nsCMatrix"
stopifnot(isSymmetric(xlx))
image(xlx, main=paste("crossprod(nmm) : Sparse", class(xlx)))
}
\keyword{classes}
\keyword{algebra}
