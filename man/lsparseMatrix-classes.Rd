\name{lsparseMatrix-classes}
\docType{class}
\alias{lsparseMatrix-class}
\alias{lgCMatrix-class}
\alias{ltCMatrix-class}
\alias{lsCMatrix-class}
\alias{lgRMatrix-class}
\alias{ltRMatrix-class}
\alias{lsRMatrix-class}
\alias{lgTMatrix-class}
\alias{ltTMatrix-class}
\alias{lsTMatrix-class}
%
\alias{-,lsparseMatrix,missing-method}
\alias{!,lsparseMatrix-method}
\alias{\%*\%,lgCMatrix,lgCMatrix-method}
\alias{\%*\%,lsparseMatrix,lsparseMatrix-method}
\alias{\%*\%,lsparseMatrix,ldenseMatrix-method}
\alias{\%*\%,ldenseMatrix,lsparseMatrix-method}
\alias{chol,lsCMatrix,logical-method}
\alias{chol,lsCMatrix,missing-method}

\alias{coerce,lgCMatrix,dgCMatrix-method}
\alias{coerce,lgCMatrix,lgTMatrix-method}
\alias{coerce,lgCMatrix,lgeMatrix-method}
\alias{coerce,lgCMatrix,matrix-method}
\alias{coerce,lgTMatrix,dgTMatrix-method}
\alias{coerce,lgTMatrix,lgCMatrix-method}
\alias{coerce,lgTMatrix,matrix-method}
\alias{coerce,lsCMatrix,dgTMatrix-method}
\alias{coerce,lsCMatrix,dsCMatrix-method}
\alias{coerce,lsCMatrix,lgCMatrix-method}
\alias{coerce,lsCMatrix,lsTMatrix-method}
\alias{coerce,lsCMatrix,matrix-method}
\alias{coerce,ltCMatrix,dtCMatrix-method}
\alias{coerce,ltCMatrix,dMatrix-method}
\alias{coerce,ltCMatrix,lgCMatrix-method}
\alias{coerce,ltCMatrix,matrix-method}
\alias{coerce,matrix,lgCMatrix-method}
\alias{coerce,matrix,lgTMatrix-method}
\alias{coerce,matrix,ltCMatrix-method}
%
\alias{coerce,ltTMatrix,dtTMatrix-method}
\alias{coerce,ltTMatrix,lgTMatrix-method}
\alias{coerce,ltTMatrix,matrix-method}
%
\alias{crossprod,lgCMatrix,missing-method}
\alias{crossprod,lgTMatrix,missing-method}
\alias{crossprod,lsparseMatrix-method}
\alias{crossprod,lsparseMatrix,missing-method}
\alias{crossprod,lsparseMatrix,lsparseMatrix-method}
\alias{crossprod,lsparseMatrix,ldenseMatrix-method}
\alias{crossprod,ldenseMatrix,lsparseMatrix-method}
\alias{diag,lgCMatrix-method}
\alias{diag,lsparseMatrix-method}
\alias{image,lgCMatrix-method}
\alias{image,lsCMatrix-method}
\alias{image,ltCMatrix-method}
\alias{image,ltTMatrix-method}
\alias{t,lgCMatrix-method}
\alias{t,lgTMatrix-method}
\alias{t,lsCMatrix-method}
\alias{t,ltCMatrix-method}
\alias{tcrossprod,lgCMatrix,missing-method}
\alias{tcrossprod,lgTMatrix,missing-method}
\alias{tcrossprod,lsparseMatrix,missing-method}
\title{Sparse logical matrices}
\description{The \code{lsparseMatrix} class is a virtual class of sparse
  matrices with \code{TRUE}/\code{FALSE} entries.  Only the positions of the
  elements that are \code{TRUE} are stored.  These can be stored in the
  ``triplet'' form (classes \code{lgTMatrix}, \code{lsTMatrix}, and
  \code{ltTMatrix} which really contain pairs, not triplets) or in
  compressed column-oriented form (classes \code{lgCMatrix},
  \code{lsCMatrix}, and \code{ltCMatrix}) or in compressed row-oriented
  form (classes \code{lgRMatrix}, \code{lsRMatrix}, and
  \code{ltRMatrix}).  The second letter in the name of these non-virtual
  classes indicates \code{g}eneral, \code{s}ymmetric, or \code{t}riangular.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("lgCMatrix",
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
    \item{coerce}{\code{signature(from = "dgCMatrix", to = "lgCMatrix")}}
    \item{t}{\code{signature(x = "lgCMatrix")}: returns the transpose
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
m <- Matrix(c(0,0,2:0), 3,5, dimnames=list(LETTERS[1:3],NULL))
(dm <- as(m, "dgCMatrix"))# no dimnames for sparse here
(lm <- as(dm, "lgCMatrix"))
str(lm) # no 'x' slot
!lm     # no longer sparse

data(KNex)
lmm <- as(KNex $ mm, "lgCMatrix")
xlx <- crossprod(lmm)
image(xlx, main=paste("crossprod(lmm) : Sparse", class(xlx)))
}
\keyword{classes}
\keyword{algebra}
