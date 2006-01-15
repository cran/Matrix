\name{dgTMatrix-class}
\docType{class}
\alias{dgTMatrix-class}
\alias{+,dgTMatrix,dgTMatrix-method}
\alias{coerce,dgTMatrix,dgCMatrix-method}
\alias{coerce,dgTMatrix,dgeMatrix-method}
\alias{coerce,dgTMatrix,dsCMatrix-method}
\alias{coerce,dgTMatrix,dtTMatrix-method}
\alias{coerce,dgTMatrix,matrix-method}
\alias{coerce,dgeMatrix,dgTMatrix-method}
\alias{coerce,matrix,dgTMatrix-method}
\alias{crossprod,dgTMatrix,missing-method}
\alias{crossprod,dgTMatrix,matrix-method}
\alias{crossprod,dgTMatrix,numeric-method}
\alias{image,dgTMatrix-method}
\alias{kronecker,dgTMatrix,dgTMatrix-method}
\alias{t,dgTMatrix-method}
\title{Sparse matrices in triplet form}
\description{The \code{"dgTMatrix"} class is the class of sparse
  matrices stored as (possibly redundant) triplets.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("dgTMatrix", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{i}:}{Object of class \code{"integer"} - the row indices
      of non-zero entries.}
    \item{\code{j}:}{Object of class \code{"integer"} - the column
      indices of non-zero entries.  Must be the same length as slot \code{i}.}
    \item{\code{x}:}{Object of class \code{"numeric"} - the (non-zero)
      entry at position \code{(i,j)}.  Must be the same length as slot
      \code{i}. If an index pair occurs more than once the corresponding
      values of slot \code{x} are added to form the element of the matrix.}
    \item{\code{Dim}:}{Object of class \code{"integer"} of length 2 -
      the dimensions of the matrix.}
  }
}
\section{Methods}{
  \describe{
    \item{+}{\code{signature(e1 = "dgTMatrix", e2 = "dgTMatrix")}}
    \item{coerce}{\code{signature(from = "dgTMatrix", to = "dgCMatrix")}}
    \item{coerce}{\code{signature(from = "dgTMatrix", to = "dgeMatrix")}}
    \item{coerce}{\code{signature(from = "dgTMatrix", to = "matrix")}}
    \item{coerce}{\code{signature(from = "dgTMatrix", to = "dsCMatrix")}}
    \item{coerce}{\code{signature(from = "matrix", to = "dsTMatrix")}}
    \item{crossprod}{\code{signature(x = "dgTMatrix", y = "missing")}
      returns \code{t(x) \%*\% x} as an \code{dsCMatrix} object.}
    \item{crossprod}{\code{signature(x = "dgTMatrix", y = "matrix")}: ... }
    \item{crossprod}{\code{signature(x = "dgTMatrix", y = "numeric")}: ... }
    \item{image}{\code{signature(x = "dgTMatrix")}: plots an image of
      \code{x} using the \code{\link[lattice]{levelplot}} function}
    \item{t}{\code{signature(x = "dgTMatrix")}: returns the transpose of
      \code{x}}
  }
}
%\references{}
%\author{}
\note{Triplet matrices are a convenient form in which to construct sparse
  matrices after which they can be coerced to \code{dgCMatrix} objects.}
\seealso{
  \code{\link{dgCMatrix-class}}
}
\examples{
m <- Matrix(1:28, nrow = 4)
##FIXME: m[c()...] <- 0
(mT <- as(m, "dgTMatrix"))

mT[1,]
mT[4, drop = FALSE]
}
\keyword{classes}
\keyword{algebra}
