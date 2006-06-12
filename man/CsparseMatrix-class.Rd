\name{CsparseMatrix-class}
\docType{class}
\alias{CsparseMatrix-class}
\alias{\%*\%,CsparseMatrix,CsparseMatrix-method}
\alias{\%*\%,CsparseMatrix,denseMatrix-method}
\alias{coerce,CsparseMatrix,TsparseMatrix-method}
\alias{coerce,CsparseMatrix,denseMatrix-method}
\alias{crossprod,CsparseMatrix,missing-method}
\alias{t,CsparseMatrix-method}
\alias{tcrossprod,CsparseMatrix,missing-method}
\alias{colMeans,CsparseMatrix-method}
\alias{colSums,CsparseMatrix-method}
\alias{rowMeans,CsparseMatrix-method}
\alias{rowSums,CsparseMatrix-method}
\title{Class "CsparseMatrix" of Sparse Matrices in Column-compressed Form}
\description{The \code{"CsparseMatrix"} class is the virtual class of
  all sparse matrices coded in sorted compressed column-oriented form.
  Since it is a virtual class, no objects may be created from it.  See
  \code{showClass("CsparseMatrix")} for its subclasses.
}

\section{Slots}{
  \describe{
    \item{\code{i}:}{Object of class \code{"integer"} of length nnzero
      (number of non-zero elements).  These are the row numbers for
      each non-zero element in the matrix.}
    \item{\code{p}:}{Object of class \code{"integer"} of pointers, one
      for each column, to the initial (zero-based) index of elements in
      the column.}
    \item{\code{factors}, \code{Dim}, \code{Dimnames}:}{inherited from
      the superclass, see the \code{\linkS4class{sparseMatrix}} class.}
  }
}
\section{Extends}{
Class \code{"sparseMatrix"}, directly.
Class \code{"Matrix"}, by class \code{"sparseMatrix"}.
}
\section{Methods}{
  \describe{
    \item{crossprod}{\code{signature(x = "CsparseMatrix", y = "missing")}: ... }
    \item{t}{\code{signature(x = "CsparseMatrix")}: ... }
    \item{tcrossprod}{\code{signature(x = "CsparseMatrix", y = "missing")}: ... }
    \item{\%*\%}{\code{signature(x = "CsparseMatrix", y =
	"CsparseMatrix")}: ...}
    \item{\%*\%}{\code{signature(x = "CsparseMatrix", y = "denseMatrix")}: ...}
  }
}
\seealso{
  its superclass, \code{\linkS4class{sparseMatrix}}, and, e.g.,
  class \code{\linkS4class{dgCMatrix}} for the links to other classes.
}
\examples{
showClass("CsparseMatrix")
}
\keyword{classes}
