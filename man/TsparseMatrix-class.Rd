\name{TsparseMatrix-class}
\title{Class "TsparseMatrix" of Sparse Matrices in Triplet Form}
\docType{class}
\alias{TsparseMatrix-class}
%% Methods:
\alias{coerce,TsparseMatrix,CsparseMatrix-method}
% "[" are in ./Xtrct-methods.Rd
\alias{crossprod,TsparseMatrix,missing-method}
\alias{tcrossprod,TsparseMatrix,missing-method}
%
\description{The \code{"TsparseMatrix"} class is the virtual class of
  all sparse matrices coded in triplet form.  Since it is a virtual class,
  no objects may be created from it.  See
  \code{showClass("TsparseMatrix")} for its subclasses.
}
\section{Slots}{
  \describe{
    \item{\code{Dim}, \code{Dimnames}:}{from the \code{"Matrix"} class,}
    \item{\code{factors}:}{from class \code{"sparseMatrix"}, see
      \code{\link{sparseMatrix-class}}.}
    \item{\code{i}:}{Object of class \code{"integer"} - the row indices
      of non-zero entries.}
    \item{\code{j}:}{Object of class \code{"integer"} - the column
      indices of non-zero entries.  Must be the same length as slot \code{i}.}
  }
}
\section{Extends}{
Class \code{"sparseMatrix"}, directly.
Class \code{"Matrix"}, by class \code{"sparseMatrix"}.
}
\section{Methods}{
  Extraction (\code{"["}) methods, see
  \code{\link{[-methods}}.%-> ./Xtrct-methods.Rd
}
% \author{Martin Maechler}
\seealso{
  its superclass, \code{\link{sparseMatrix-class}}, and
  \code{\link{dgTMatrix-class}} for the  links to other classes
}
\examples{
showClass("TsparseMatrix")
## or just the subclasses' names
names(getClass("TsparseMatrix")@subclasses)
}
\keyword{classes}
