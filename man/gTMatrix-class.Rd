\name{gTMatrix-class}
\title{Class "gTMatrix" of General Sparse Matrices in Triplet Form}
\docType{class}
\alias{gTMatrix-class}
% "[" are in ./Xtrct-methods.Rd
%
\description{  ~~ A concise (1-5 lines) description of what the class is.  ~~}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{Dim}:}{class \code{"integer"} of length 2, as for
      \code{\link{Matrix-class}}.}
    \item{\code{Dimnames}:}{class \code{"list"} as in the \code{"Matrix"} class.}
  }
}
\section{Extends}{
Class \code{"sparseMatrix"}, directly.
Class \code{"Matrix"}, by class \code{"sparseMatrix"}.
}
\section{Methods}{
  Extraction (\code{"["}) methods, see \code{\link{[-methods}}.
}
% \author{Martin Maechler}
\seealso{
  its superclass, \code{\link{sparseMatrix-class}}, and
  \code{\link{dgTMatrix-class}} for the  links to other classes
}
% \examples{
% }
\keyword{classes}
