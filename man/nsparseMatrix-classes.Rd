\name{nsparseMatrix-classes}
\docType{class}
\alias{nsparseMatrix-class}
\alias{ngCMatrix-class}
\alias{ntCMatrix-class}
\alias{nsCMatrix-class}
\alias{ngRMatrix-class}
\alias{ntRMatrix-class}
\alias{nsRMatrix-class}
\alias{ngTMatrix-class}
\alias{ntTMatrix-class}
\alias{nsTMatrix-class}
%
\alias{-,nsparseMatrix,missing-method}
\alias{!,nsparseMatrix-method}
\alias{\%*\%,ngCMatrix,ngCMatrix-method}
\alias{\%*\%,nsparseMatrix,nsparseMatrix-method}
\alias{\%*\%,nsparseMatrix,ndenseMatrix-method}
\alias{\%*\%,ndenseMatrix,nsparseMatrix-method}
\alias{chol,nsCMatrix,logical-method}
\alias{chol,nsCMatrix,missing-method}

\alias{coerce,ngCMatrix,dgCMatrix-method}
\alias{coerce,ngCMatrix,ngTMatrix-method}
\alias{coerce,ngCMatrix,ngeMatrix-method}
\alias{coerce,ngCMatrix,matrix-method}
\alias{coerce,ngTMatrix,dgTMatrix-method}
\alias{coerce,ngTMatrix,ngCMatrix-method}
\alias{coerce,ngTMatrix,matrix-method}
\alias{coerce,nsCMatrix,dgTMatrix-method}
\alias{coerce,nsCMatrix,dsCMatrix-method}
\alias{coerce,nsCMatrix,ngCMatrix-method}
\alias{coerce,nsCMatrix,nsTMatrix-method}
\alias{coerce,nsCMatrix,matrix-method}
\alias{coerce,nsTMatrix,dsTMatrix-method}
\alias{coerce,nsTMatrix,ngTMatrix-method}
\alias{coerce,nsTMatrix,matrix-method}
\alias{coerce,ntCMatrix,dtCMatrix-method}
\alias{coerce,ntCMatrix,dMatrix-method}
\alias{coerce,ntCMatrix,ngCMatrix-method}
\alias{coerce,ntCMatrix,matrix-method}
\alias{coerce,matrix,ngCMatrix-method}
\alias{coerce,matrix,ngTMatrix-method}
\alias{coerce,matrix,ntCMatrix-method}
%
\alias{coerce,ntTMatrix,dtTMatrix-method}
\alias{coerce,ntTMatrix,ngTMatrix-method}
\alias{coerce,ntTMatrix,matrix-method}
\alias{coerce,matrix,ntTMatrix-method}
%
\alias{crossprod,ngCMatrix,missing-method}
\alias{crossprod,ngTMatrix,missing-method}
\alias{crossprod,nsparseMatrix-method}
\alias{crossprod,nsparseMatrix,missing-method}
\alias{crossprod,nsparseMatrix,nsparseMatrix-method}
\alias{crossprod,nsparseMatrix,ndenseMatrix-method}
\alias{crossprod,ndenseMatrix,nsparseMatrix-method}
\alias{diag,ngCMatrix-method}
\alias{diag,nsparseMatrix-method}
\alias{image,ngCMatrix-method}
\alias{image,nsCMatrix-method}
\alias{image,nsTMatrix-method}
\alias{image,ntCMatrix-method}
\alias{image,ntTMatrix-method}
\alias{t,ngCMatrix-method}
\alias{t,ngTMatrix-method}
\alias{t,nsCMatrix-method}
\alias{t,ntCMatrix-method}
\alias{tcrossprod,ngCMatrix,missing-method}
\alias{tcrossprod,ngTMatrix,missing-method}
\alias{tcrossprod,nsparseMatrix,missing-method}
\title{Sparse logical matrices}
\description{The \code{nsparseMatrix} class is a virtual class of sparse
  matrices with \code{TRUE}/\code{FALSE} entries.  Only the positions of the
  elements that are \code{TRUE} are stored.  These can be stored in the
  ``triplet'' form (classes \code{ngTMatrix}, \code{nsTMatrix}, and
  \code{ntTMatrix} which really contain pairs, not triplets) or in
  compressed column-oriented form (classes \code{ngCMatrix},
  \code{nsCMatrix}, and \code{ntCMatrix}) or in compressed row-oriented
  form (classes \code{ngRMatrix}, \code{nsRMatrix}, and
  \code{ntRMatrix}).  The second letter in the name of these non-virtual
  classes indicates \code{g}eneral, \code{s}ymmetric, or \code{t}riangular.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("ngCMatrix",
    ...)} and so on.  More frequently objects are created by coercion of
  a numeric sparse matrix to the logical form for use in
  the symbolic analysis phase
  of an algorithm involving sparse matrices.  Such algorithms often
  involve two phases: a symbolic phase wherein the positions of the
  non-zeros in the result are determined and a numeric phase wherein the
  actual results are calculated.  During the symbolic phase only the
  positions of the non-zero elements in any operands are of interest,
  hence any numeric sparse matrices can be treated as logical sparse
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
      Present in triplet and compressed column-oriented forms only.}
    \item{\code{Dim}:}{Object of class \code{"integer"} - the dimensions
      of the matrix.}
  }
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "dgCMatrix", to = "ngCMatrix")}}
    \item{t}{\code{signature(x = "ngCMatrix")}: returns the transpose
      of \code{x}}
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
nm <- as(m, "ngCMatrix")
str(nm) # no 'x' slot
!nm     # no longer sparse
## consistency check:
stopifnot(xor(as( nm, "matrix"),
              as(!nm, "matrix")))

data(KNex)
nmm <- as(KNex $ mm, "ngCMatrix")
str(xlx <- crossprod(nmm))# "nsCMatrix"
stopifnot(isSymmetric(xlx))
image(xlx, main=paste("crossprod(nmm) : Sparse", class(xlx)))
}
\keyword{classes}
\keyword{algebra}
