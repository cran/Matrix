\name{sscCrosstab-class}
\docType{class}
\alias{sscCrosstab-class}

\title{Pairwise crosstabulations}
\description{An \code{sscCrosstab} object is a pairwise crosstabulation
  of two or more factors expressed as an \code{sscMatrix} plus
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
    \item{\code{uplo}:}{as for \code{"sscMatrix"} objects.}
    \item{\code{p}:}{as for \code{"sscMatrix"} objects.}
    \item{\code{i}:}{as for \code{"sscMatrix"} objects.}
    \item{\code{x}:}{as for \code{"sscMatrix"} objects.}
    \item{\code{Dim}:}{as for \code{"sscMatrix"} objects.}
    \item{\code{factorization}:}{as for \code{"sscMatrix"} objects.}
  }
}
\section{Extends}{
Class \code{"sscMatrix"}, directly.
Class \code{"cscMatrix"}, by class \code{"sscMatrix"}.
}
%\section{Methods}{}
%\references{}
%\author{}
%\note{}
\seealso{
 \code{\link{sscMatrix-class}}, \code{\link{sscCrosstab}}
}
%\examples{}
\keyword{classes}
