\name{pdfactor-class}
\docType{class}
\alias{pdfactor-class}
\title{Factors of positive-definite matrices}
\description{If \eqn{\Sigma}{S} is a positive-definite matrix then any
  square matrix \eqn{L}{L} such that \eqn{\Sigma = L'L}{S = L'L} is said
  to be a square-root factor of \eqn{\Sigma}{S}.  The \code{pdfactor}
  class is the class of such matrices.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("pdfactor",
    ...)}.  More frequently they are created by coercing a \code{pdMat}
  or \code{pdmatrix} object to the \code{pdfactor} class.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{\code{"matrix"}: the factor itself, a numeric,
      square matrix}
    \item{\code{logDet}:}{\code{"numeric"}: the logarithm of the
      absolute value of the determinant of the factor.}
  }
}
\section{Extends}{
Class \code{"matrix"}, from data part.
Class \code{"structure"}, by class "matrix".
Class \code{"array"}, by class "matrix".
Class \code{"vector"}, by class "matrix".
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "pdIdent", to = "pdfactor")}:
      extract the pdfactor from a \code{pdIdent} object.}
  }
}
%\references{}
%\author{Douglas Bates \email{bates@stat.wisc.edu} and Saikat DebRoy \email{saikat@stat.wisc.edu}}
%\note{}
\seealso{\code{\link{pdMat-class}}, \code{\link{pdmatrix-class}}}
%\examples{}
\keyword{classes}
