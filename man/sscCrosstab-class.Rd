\name{sscCrosstab-class}
\docType{class}
\alias{sscCrosstab-class}
\alias{image,sscCrosstab-method}
\title{Pairwise crosstabulations}
\description{An \code{sscCrosstab} object is a pairwise crosstabulation
  of two or more factors expressed as an \code{dsCMatrix} plus
  additional information.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("sscCrosstab", ...)}
or by a call to the constructor function \code{sscCrosstab}.
}
\section{Slots}{
  \describe{
    \item{\code{Gp}:}{integer vector of pointers to the beginning of
      each group of columns.  Its length is 1 greater than the number of
      groups.  Successive differences give the number of levels of each
      factor.} 
    \item{\code{perm}:}{Integer permutation vector giving a
      fill-reducing permutation for the Cholesky decomposition.  The
      permutation is constrained to keep the groups of columns
      corresponding to different factors distinct.}
    \item{\code{uplo}:}{as for \code{"dsCMatrix"} objects.}
    \item{\code{p}:}{as for \code{"dsCMatrix"} objects.}
    \item{\code{i}:}{as for \code{"dsCMatrix"} objects.}
    \item{\code{x}:}{as for \code{"dsCMatrix"} objects.}
    \item{\code{Dim}:}{as for \code{"dsCMatrix"} objects.}
    \item{\code{factors}:}{as for \code{"dsCMatrix"} objects.}
  }
}
\section{Extends}{
Class \code{"dsCMatrix"}, directly.
Class \code{"dgCMatrix"}, by class \code{"dsCMatrix"}.
}
\section{Methods}{
  \describe{
    \item{image}{\code{signature(x = "sscCrosstab")}: plots an image of
      \code{x} using the \code{\link[lattice]{levelplot}} function}
  }
}
%\references{}
%\author{}
%\note{}
\seealso{
 \code{\link{dsCMatrix-class}}, \code{\link{sscCrosstab}}
}
%\examples{}
\keyword{classes}
