\name{tripletMatrix-class}
\docType{class}
\alias{tripletMatrix-class}
\alias{+,tripletMatrix,tripletMatrix-method}
\alias{coerce,tripletMatrix,geMatrix-method}
\alias{coerce,tripletMatrix,matrix-method}
\alias{coerce,tripletMatrix,cscMatrix-method}
\alias{coerce,tripletMatrix,sscMatrix-method}
\alias{image,tripletMatrix-method}
\alias{t,tripletMatrix-method}
\title{Class "tripletMatrix" sparse matrices in triplet form}
\description{The \code{"tripletMatrix"} class is the class of sparse
  matrices stored as (possibly redundant) triplets.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("tripletMatrix", ...)}.
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
    \item{+}{\code{signature(e1 = "tripletMatrix", e2 = "tripletMatrix")}}
    \item{coerce}{\code{signature(from = "tripletMatrix", to = "cscMatrix")}}
    \item{coerce}{\code{signature(from = "tripletMatrix", to = "geMatrix")}}
    \item{coerce}{\code{signature(from = "tripletMatrix", to = "matrix")}}
    \item{coerce}{\code{signature(from = "tripletMatrix", to = "sscMatrix")}}
    \item{image}{\code{signature(x = "tripletMatrix")}: plots an image of
      \code{x} using the \code{\link[lattice]{levelplot}} function}
    \item{t}{\code{signature(x = "tripletMatrix")}: returns the transpose of
      \code{x}}
  }
}
%\references{}
%\author{}
\note{Triplet matrices are a convenient form in which to construct sparse
  matrices after which they can be coerced to \code{cscMatrix} objects.}
\seealso{
  \code{\link{cscMatrix-class}}
}
%\examples{}
\keyword{classes}
\keyword{algebra}
