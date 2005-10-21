\name{CsparseMatrix-class}
\docType{class}
\alias{CsparseMatrix-class}
\alias{crossprod,CsparseMatrix,missing-method}
\alias{t,CsparseMatrix-method}
\alias{tcrossprod,CsparseMatrix-method}
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
      the superclass, see \code{\link{sparseMatrix-class}}.}
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
    \item{tcrossprod}{\code{signature(x = "CsparseMatrix")}: ... }
  }
}
\seealso{
  its superclass, \code{\link{sparseMatrix-class}}, and, e.g.,
  \code{\link{dgCMatrix-class}} for the links to other classes.
}
\examples{
showClass("CsparseMatrix")
}
\keyword{classes}
