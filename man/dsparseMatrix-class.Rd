\name{dsparseMatrix-class}
\title{Virtual Class "dsparseMatrix" of Numeric Sparse Matrices}
\docType{class}
%
\alias{dsparseMatrix-class}
%
\alias{Arith,dsparseMatrix,logical-method}
\alias{Arith,dsparseMatrix,numeric-method}
\alias{Arith,logical,dsparseMatrix-method}
\alias{Arith,numeric,dsparseMatrix-method}
\alias{Ops,dsparseMatrix,nsparseMatrix-method}
\alias{Summary,dsparseMatrix-method}
\alias{as.logical,dsparseMatrix-method}
\alias{as.numeric,dsparseMatrix-method}
\alias{coerce,dsparseMatrix,lMatrix-method}
\alias{coerce,dsparseMatrix,lsparseMatrix-method}
\alias{coerce,dsparseMatrix,nMatrix-method}
\alias{coerce,dsparseMatrix,nsparseMatrix-method}
\alias{coerce,matrix,dsparseMatrix-method}
\alias{coerce,numLike,dsparseMatrix-method}
%
\description{The Class \code{"dsparseMatrix"} is the virtual (super) class of
  all numeric sparse matrices.
}
\section{Slots}{
  \describe{
    \item{\code{Dim}:}{the matrix dimension, see class \code{"\linkS4class{Matrix}"}.}
    \item{\code{Dimnames}:}{see the \code{"Matrix"} class.}
    \item{\code{x}:}{a \code{\link{numeric}} vector containing the
      (non-zero) matrix entries.}
  }
}
\section{Extends}{
  Class \code{"dMatrix"} and \code{"sparseMatrix"}, directly.\cr
  Class \code{"Matrix"}, by the above classes.
}
% \section{Methods}{
% No methods defined with class "dsparseMatrix" in the signature.
% }
%%\author{Martin}
\seealso{
  the documentation of the (non virtual) sub classes, see
  \code{showClass("dsparseMatrix")}; in particular,
  \linkS4class{dgTMatrix}, \linkS4class{dgCMatrix}, and
  \linkS4class{dgRMatrix}.
}
\examples{
showClass("dsparseMatrix")
}
\keyword{classes}
