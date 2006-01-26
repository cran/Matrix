\name{lgeMatrix-class}
\docType{class}
\alias{lgeMatrix-class}
\alias{as.vector,lgeMatrix,missing-method}
\alias{coerce,matrix,lgeMatrix-method}
\alias{coerce,lgeMatrix,dgeMatrix-method}
\alias{coerce,lgeMatrix,matrix-method}
\alias{coerce,lgeMatrix,lgTMatrix-method}
\alias{coerce,lgeMatrix,lsyMatrix-method}
\alias{t,lgeMatrix-method}

\title{Class "lgeMatrix" of General Dense Logical Matrices}
\description{This is the class of general dense \code{\link{logical}}
  matrices.
}
\section{Slots}{
  \describe{
    \item{\code{x}:}{Object of class \code{"logical"}. The logical
      values that constitute the matrix, stored in column-major order.}
    \item{\code{Dim},\code{Dimnames}:}{The dimension (a length-2
      \code{"integer"}) and corresponding names (or \code{NULL}), see the
      \code{\link{Matrix-class}}.}
    \item{\code{factors}:}{Object of class \code{"list"}.  A named
      list of factorizations that have been computed for the matrix.}
  }
}
\section{Extends}{
Class \code{"ldenseMatrix"}, directly.
Class \code{"lMatrix"}, by class \code{"ldenseMatrix"}.
Class \code{"denseMatrix"}, by class \code{"ldenseMatrix"}.
Class \code{"Matrix"}, by class \code{"ldenseMatrix"}.
Class \code{"Matrix"}, by class \code{"ldenseMatrix"}.
}
\section{Methods}{
  Currently, mainly \code{\link{t}()} and coercion methods (for
  \code{\link{as}(.)}); use, e.g.,
  \code{\link{showMethods}(class="lgeMatrix")} for details.
}
\seealso{
  Non-general logical dense classes such as
  \code{\link{ltrMatrix-class}}, or \code{\link{lsyMatrix-class}};
  \emph{sparse} logical matrices such as \code{\link{lgCMatrix-class}}.
}
\examples{
showClass("lgeMatrix")
str(new("lgeMatrix"))
}
\keyword{classes}
