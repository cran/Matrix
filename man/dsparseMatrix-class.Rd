\name{dsparseMatrix-class}
\docType{class}
\alias{dsparseMatrix-class}
\title{Virtual Class "dsparseMatrix" of Numeric Sparse Matrices}
\description{The Class \code{"dsparseMatrix"} is the virtual (super) class of
  all numeric sparse matrices.
}
\section{Slots}{
  \describe{
    \item{\code{Dim}:}{the matrix dimension, see \link{Matrix-class}.}
    \item{\code{Dimnames}:}{see the "Matrix" class.}

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
  \link{dgTMatrix-class}, \link{dgCMatrix-class}, and
  \link{dgRMatrix-class}.
}
\examples{
showClass("dsparseMatrix")
}
\keyword{classes}
