\name{pdmatrix-class}
\docType{class}
\alias{pdmatrix-class}
\alias{coerce,matrix,pdmatrix-method}
\title{Positive-definite matrices }
\description{The class of positive-definite, symmetric matrices.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("pdmatrix", ...)}.
  More frequently they are created by coercing a \code{pdMat}
  object to the \code{pdmatrix} class.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"matrix"}, which should
      be a numeric matrix.}
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
%     \item{coerce}{\code{signature(from = "pdDiag", to = "pdmatrix")}: extract the positive-definite matrix }
%     \item{coerce}{\code{signature(from = "pdIdent", to = "pdmatrix")}:
%       extract the positive-definite matrix}
%     \item{coerce}{\code{signature(from = "pdLogChol", to = "pdmatrix")}:
%       extract the positive-definite matrix}
%     \item{coerce}{\code{signature(from = "pdNatural", to = "pdmatrix")}:
%       extract the positive-definite matrix}
%    \item{coerce}{\code{signature(from = "pdScalar", to = "pdmatrix")}:
%      extract the positive-definite matrix}
    \item{coerce}{\code{signature(from = "matrix", to = "pdmatrix")}:
      create a \code{pdmatrix} object from a matrix (does not check for
      symmetry or positive-definiteness)}
  }
}
%\references{}
%\author{Douglas Bates \email{bates@stat.wisc.edu} and Saikat DebRoy \email{saikat@stat.wisc.edu}}
%\note{}
\seealso{\code{\link{dpoMatrix-class}}}}
%\examples{}
\keyword{classes}
