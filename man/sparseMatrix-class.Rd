\name{sparseMatrix-class}
\docType{class}
\alias{sparseMatrix-class}
\alias{show,sparseMatrix-method}
\alias{\%*\%,sparseMatrix,ddenseMatrix-method}
\alias{\%*\%,ddenseMatrix,sparseMatrix-method}
\alias{crossprod,sparseMatrix,ddenseMatrix-method}
\alias{crossprod,ddenseMatrix,sparseMatrix-method}
\alias{coerce,graphNEL,sparseMatrix-method}
\title{Virtual Class "sparseMatrix" --- Mother of Sparse Matrices}
\description{Virtual Mother Class of All Sparse Matrices}
\section{Slots}{
  \code{Dim}, and \code{Dimnames}, see \code{\link{Matrix-class}}.
% Should this be "written in stone" or should they be added?
%   Note that the current implementation keeps \code{Dimnames} empty for
%   sparse matrices.
}
\section{Extends}{
  Class \code{"Matrix"}, directly.
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "sparseMatrix")}: The
      \code{\link{show}} method for sparse matrices prints
      \emph{\dQuote{structural}} zeroes as \code{"."} or blank.}
  }
}
\note{
  In method selection for multiplication operations (i.e. \code{\%*\%}
  and the two-argument form of \code{\link[base]{crossprod}})
  the sparseMatrix class takes precedence in the sense that if one
  operand is a sparse matrix and the other is any type of dense matrix
  then the dense matrix is coerced to a \code{dgeMatrix} and the
  appropriate sparse matrix method is used.
}
%\author{Martin}
\examples{
showClass("sparseMatrix") ## and look at the help() of its subclasses
}
\keyword{classes}
