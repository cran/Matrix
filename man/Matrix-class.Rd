\name{Matrix-class}
\docType{class}
\alias{Matrix-class}
\alias{show,Matrix-method}
\title{Class "Matrix" Virtual Class of Matrices}
\description{
  The \code{Matrix} class is a virtual class contained by all actual
  classes in the \pkg{Matrix} package.
}
%% FIXME: It's no longer virtual in the strict sense
%\section{Objects from the Class}{A virtual Class: No objects may be
%  created from it.
%}
\section{Slots}{
  Common to \emph{all} matrix object in the package:
  \describe{
    \item{\code{Dim}:}{Object of class \code{"integer"} - the dimensions
     of the matrix - must be an integer vector with exactly two
     non-negative values.}
    \item{\code{Dimnames}:}{list of length two; each component
      containing NULL or a \code{\link{character}} vector length
      equal the corresponding \code{Dim} element.}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "Matrix"}}
  }
}
%\references{}
\author{Douglas Bates \email{bates@stat.wisc.edu}}
\seealso{
  \code{\link{dgeMatrix-class}}, \code{\link{dgCMatrix-class}}, and
  \code{\link{Matrix}} for construction (and examples).
}
\examples{
 (M <- new("Matrix"))
 (M <- new("Matrix", Dim = 1:2, Dimnames = list("r1", c("c1","c2"))))
 str(M)
}
\keyword{classes}
\keyword{algebra}
