\name{RsparseMatrix-class}
\docType{class}
\alias{RsparseMatrix-class}
%% No methods currently !
\title{Class "RsparseMatrix" of Sparse Matrices in Column-compressed Form}
\description{The \code{"RsparseMatrix"} class is the virtual class of
  all sparse matrices coded in sorted compressed row-oriented form.
  Since it is a virtual class, no objects may be created from it.  See
  \code{showClass("RsparseMatrix")} for its subclasses.
}
\section{Slots}{
  \describe{
    \item{\code{j}:}{Object of class \code{"integer"} of length \code{nnzero}
      (number of non-zero elements).  These are the row numbers for
      each non-zero element in the matrix.}
    \item{\code{p}:}{Object of class \code{"integer"} of pointers, one
      for each row, to the initial (zero-based) index of elements in
      the row.}
    \item{\code{Dim}, \code{Dimnames}:}{inherited from
      the superclass, see \code{\link{sparseMatrix-class}}.}
  }
}
\section{Extends}{
  Class \code{"sparseMatrix"}, directly.
  Class \code{"Matrix"}, by class \code{"sparseMatrix"}.
}
\section{Methods}{
  \bold{NO} methods are defined currently, since we rather use the
  \code{\link{CsparseMatrix-class}} in \pkg{Matrix}.
%   \describe{
%     \item{crossprod}{\code{signature(x = "RsparseMatrix", y = "missing")}: ... }
%     \item{t}{\code{signature(x = "RsparseMatrix")}: ... }
%     \item{tcrossprod}{\code{signature(x = "RsparseMatrix", y = "missing")}: ... }
%     \item{\%*\%}{\code{signature(x = "RsparseMatrix", y =
% 	"RsparseMatrix")}: ...}
%     \item{\%*\%}{\code{signature(x = "RsparseMatrix", y = "denseMatrix")}: ...}
%   }
}
\seealso{
  its superclass, \code{\link{sparseMatrix-class}}, and, e.g.,
  \code{\link{dgRMatrix-class}} for the links to other classes.
}
\examples{
showClass("RsparseMatrix")
}
\keyword{classes}
