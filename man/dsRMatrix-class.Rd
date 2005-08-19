\name{dsRMatrix-class}
\title{Symmetric Sparse Compressed Row Matrices}
\docType{class}
\alias{dsRMatrix-class}
%%--- No methods currently --- see ./dsRMatrix-class.Rd for more
\description{The \code{dsRMatrix} class is a class of symmetric, sparse
  matrices in the compressed, row-oriented format.  In this
  implementation the non-zero elements in the rows are sorted into
  increasing row order.
}
\section{Objects from the Class}{
  These classes are currently still mostly unimplemented!

  Objects can be created by calls of the form \code{new("dsRMatrix", ...)}
  or \code{new("dtRMatrix", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{uplo}:}{A character object indicating if the upper
      triangle (\code{"U"} or \code{"u"}) or the lower triangle
      (\code{"L"} or \code{"l"}) is stored.  At present only the lower
      triangle form is allowed.}
%% FIXME: this was by copy-paste-edit: is it ok?
    \item{\code{j}:}{Object of class \code{"integer"} of length
      \code{nnzero} (number of non-zero elements).  These are the row
      numbers for each non-zero element in the matrix.}
    \item{\code{p}:}{Object of class \code{"integer"} of pointers, one
      for each row, to the initial (zero-based) index of elements in
      the row. (Only present in the \code{dsRMatrix} class.)}
    \item{\code{factors}:}{Object of class \code{"list"} - a list
      of factorizations of the matrix. }
    \item{\code{x}:}{Object of class \code{"numeric"} - the non-zero
      elements of the matrix.}
    \item{\code{Dim}:}{Object of class \code{"integer"} - the dimensions
     of the matrix - must be an integer vector with exactly two
     non-negative values.}
    \item{\code{Dimnames}:}{List of length two, see \code{\link{Matrix}}.}
  }
}
\section{Extends}{
Class \code{"dgRMatrix"}, directly.
Class \code{"dsparseMatrix"}, by class \code{"dgRMatrix"}.
Class \code{"dMatrix"}, by class \code{"dgRMatrix"}.
Class \code{"sparseMatrix"}, by class \code{"dgRMatrix"}.
Class \code{"Matrix"}, by class \code{"dgRMatrix"}.
Class \code{"Matrix"}, by class \code{"dgRMatrix"}.
}
\section{Methods}{
  No methods currently with class "dsRMatrix" in the signature.
%   \describe{
%     \item{solve}{\code{signature(a = "dsRMatrix", b = "matrix")}: Solve
%       a linear system of equations defined by \code{x} using a Cholesky
%       decomposition.}
% ......
%     \item{coerce}{\code{signature(from = "dsRMatrix", to = "dgTMatrix")}}
% ......
%   }
}
%\references{}
%\author{}
%\note{}
\seealso{
  \code{\link{dgCMatrix-class}}, \code{\link{dgTMatrix-class}},
  \code{\link{dgeMatrix-class}}
}
\examples{
m0 <- new("dsRMatrix")
## FIXME: print() fails
m2 <- new("dsRMatrix", x = 1:4, Dim = as.integer(c(2,2)))
## FIXME: print() fails
str(m2)
if(FALSE)## FIXME
str(m3 <- as(Matrix(diag(2)), "dsRMatrix"))
}
\keyword{classes}
\keyword{algebra}
