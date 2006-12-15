\name{dsCMatrix-class}
\title{Numeric Symmetric Sparse (column compressed) Matrices}
\docType{class}
\alias{dsCMatrix-class}
\alias{dsTMatrix-class}
%
\alias{solve,dsCMatrix,dsparseMatrix-method}
\alias{solve,dsCMatrix,matrix-method}
\alias{solve,dsCMatrix,numeric-method}
\alias{solve,dsCMatrix,ddenseMatrix-method}
\alias{chol,dsCMatrix,missing-method}
\alias{chol,dsCMatrix,logical-method}
\alias{t,dsCMatrix-method}
\alias{t,dsTMatrix-method}
\alias{coerce,dgeMatrix,dsCMatrix-method}
\alias{coerce,dsCMatrix,dgCMatrix-method}
\alias{coerce,dsCMatrix,dgTMatrix-method}
\alias{coerce,dsCMatrix,dgeMatrix-method}
\alias{coerce,dsCMatrix,dsTMatrix-method}
\alias{coerce,dsCMatrix,dsyMatrix-method}
\alias{coerce,dsCMatrix,lsCMatrix-method}
\alias{coerce,dsCMatrix,nsCMatrix-method}
\alias{coerce,dsCMatrix,matrix-method}
\alias{coerce,dgeMatrix,dsTMatrix-method}
\alias{coerce,dsTMatrix,dgTMatrix-method}
\alias{coerce,dsTMatrix,dgeMatrix-method}
\alias{coerce,dsTMatrix,dsCMatrix-method}
\alias{coerce,dsTMatrix,dsyMatrix-method}
\alias{coerce,dsTMatrix,lsTMatrix-method}
\alias{coerce,dsTMatrix,matrix-method}
\alias{coerce,dsyMatrix,dsCMatrix-method}
\alias{coerce,dsyMatrix,dsTMatrix-method}
\alias{coerce,matrix,dsCMatrix-method}
\alias{coerce,matrix,dsTMatrix-method}
\alias{coerce,graphNEL,dsCMatrix-method}
\alias{determinant,dsCMatrix,missing-method}
\alias{determinant,dsCMatrix,logical-method}
\alias{image,dsCMatrix-method}
\description{The \code{dsCMatrix} class is a class of symmetric, sparse
  numeric matrices in the compressed, column-oriented format.  In this
  implementation the non-zero elements in the columns are sorted into
  increasing row order.\cr
  The \code{dsTMatrix} class is the class of
  symmetric, sparse numeric matrices in triplet format.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("dsCMatrix",
    ...)} or \code{new("dsTMatrix", ...)}
}
\section{Slots}{
  \describe{
    \item{\code{uplo}:}{A character object indicating if the upper
      triangle (\code{"U"} or \code{"u"}) or the lower triangle
      (\code{"L"} or \code{"l"}) is stored.}
    \item{\code{p}:}{Object of class \code{"integer"} of pointers, one
      for each column, to the initial (zero-based) index of elements in
      the column. (Only present in the \code{dsCMatrix} class.)}
    \item{\code{i}:}{Object of class \code{"integer"} of length nnzero
      (half number of non-zero elements).  These are the row numbers for
      each non-zero element in the lower triangle of the matrix.}
    \item{\code{j}:}{Object of class \code{"integer"} of length nnzero
      (half number of non-zero elements).  These are the column numbers for
      each non-zero element in the lower triangle of the matrix. (Only
      present in the \code{dsTMatrix} class.)}
    \item{\code{x}:}{Object of class \code{"numeric"} - the non-zero
      elements of the matrix (to be duplicated for full matrix).}
    \item{\code{factors}:}{Object of class \code{"list"} - a list
      of factorizations of the matrix. }
    \item{\code{Dim}:}{Object of class \code{"integer"} - the dimensions
     of the matrix - must be an integer vector with exactly two
     non-negative values.}
  }
}
\section{Extends}{
Class \code{"dgCMatrix"}, directly.
}
\section{Methods}{
  \describe{
    \item{solve}{\code{signature(a = "dsCMatrix", b = "dsparseMatrix")}: Solve
      a linear system of equations defined by \code{x} using a Cholesky
      decomposition.  All steps will be based on \emph{sparse}
      representations.}
    \item{solve}{\code{signature(a = "dsCMatrix", b = "matrix")}: Solve
      a linear system of equations defined by \code{x} using a Cholesky
      decomposition.}
    \item{chol}{\code{signature(x = "dsCMatrix", pivot = "logical")}:
      Returns (and stores) the Cholesky decomposition of the matrix
      \code{x}. If \code{pivot} is \code{TRUE} (the default) Metis is
      used to create a reordering of the rows and columns of \code{x} so
      as to minimize fill-in.}
    \item{determinant}{\code{signature(x = "dsCMatrix", logarithm =
	"missing")}: Evaluate the determinant of \code{x} on the
      logarithm scale.  This creates and stores the Cholesky factorization.}
    \item{determinant}{\code{signature(x = "dsCMatrix", logarithm =
	"logical")}: Evaluate the determinant of \code{x} on the
      logarithm scale or not, according to the \code{logarithm}
      argument.  This creates and stores the Cholesky factorization.}
    \item{t}{\code{signature(x = "dsCMatrix")}: Transpose.  Because
      \code{x} is symmetric this has no effect.}
    \item{t}{\code{signature(x = "dsTMatrix")}: Transpose.  For the
      \code{dsTMatrix} class the row and column indices are interchanged
      so that a matrix for which the upper triangle is stored produces a
      matrix for which the lower triangle is stored and vice versa.}
    \item{coerce}{\code{signature(from = "dsCMatrix", to = "dgTMatrix")}}
    \item{coerce}{\code{signature(from = "dsCMatrix", to = "dgeMatrix")}}
    \item{coerce}{\code{signature(from = "dsCMatrix", to = "matrix")}}
    \item{coerce}{\code{signature(from = "dsTMatrix", to = "dgeMatrix")}}
    \item{coerce}{\code{signature(from = "dsTMatrix", to = "dsCMatrix")}}
    \item{coerce}{\code{signature(from = "dsTMatrix", to = "dsyMatrix")}}
    \item{coerce}{\code{signature(from = "dsTMatrix", to = "matrix")}}
  }
}
%\references{}
%\author{}
%\note{}
\seealso{
  Classes \code{\linkS4class{dgCMatrix}}, \code{\linkS4class{dgTMatrix}},
  \code{\linkS4class{dgeMatrix}}
}
\examples{
## First a "dgCMatrix"
mm <- Matrix(toeplitz(c(10, 0, 1, 0, 3)), sparse = TRUE)
mT <- as(mm, "dgTMatrix")
(symM <- as(mT, "dsCMatrix"))
str(symM)
sym2 <- as(symM, "TsparseMatrix")
stopifnot(class(sym2) == "dsTMatrix",
	  identical(sym2[1,], symM[1,]),
	  identical(sym2[,2], symM[,2]))
}
\keyword{classes}
\keyword{algebra}
