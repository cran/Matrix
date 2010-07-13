\name{predModule-class}
\Rdversion{1.1}
\title{Class "predModule" and SubClasses}
\docType{class}
\alias{predModule-class}
\alias{dPredModule-class}
\alias{sPredModule-class}
%
\alias{coerce,ddenseModelMatrix,predModule-method}
\alias{coerce,dsparseModelMatrix,predModule-method}
\description{
  The class \code{"predModule"} and notably its subclasses
  \code{"dPredModule"} and \code{"sPredModule"} encapsulate information
  about linear predictors in statistical models.  They incorporate a
  \code{\linkS4class{modelMatrix}}, the corresponding coefficients and a
  representation of a triangular factor from the, possibly weighted or
  otherwise modified, model matrix.
}
\section{Objects from the Classes}{
  Objects are typically created by coercion from objects of class
  \code{\linkS4class{ddenseModelMatrix}} or
  \code{\linkS4class{dsparseModelMatrix}}.
}
\section{Slots}{
  The virtual class \code{"predModule"} and its two subclasses all have slots
  \describe{
    \item{\code{X}:}{a \code{\linkS4class{modelMatrix}}.}
    \item{\code{coef}:}{\code{"numeric"} coefficient vector of length
      \code{ncol(.)}\eqn{:= p}.}
    \item{\code{Vtr}:}{\code{"numeric"} vector of length \eqn{p},
      to contain \eqn{V'r} (\dQuote{\bold{V} \bold{t}ransposed \bold{r}}).}
    \item{\code{fac}:}{a representation of a triangular factor, the
      Cholesky decomposition of \eqn{V'V}.}% << FIXME? (weights !)
  }
  The actual classes \code{"dPredModule"} and \code{"sPredModule"}
  specify specific (sub) classes for the two non-trivial slots,
  \describe{
    \item{\code{X}:}{a \code{"\linkS4class{ddenseModelMatrix}"} or
      \code{"\linkS4class{dsparseModelMatrix}"}, respectively.}
    \item{\code{fac}:}{For the \code{"dpredModule"} class this factor is a
      \code{\linkS4class{Cholesky}} object.  For the \code{"spredModule"} class
      it is of class \code{\linkS4class{CHMfactor}}.}
  }
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "ddenseModelMatrix", to = "predModule")}:
      Creates a \code{"dPredModule"} object.}
    \item{coerce}{\code{signature(from = "dsparseModelMatrix", to = "predModule")}:
      Creates an \code{"sPredModule"} object.}
  }
}
\author{Douglas Bates}
\seealso{
  \code{\link{model.Matrix}()} which returns a
  \code{"\linkS4class{ddenseModelMatrix}"} or
  \code{"\linkS4class{dsparseModelMatrix}"} object, depending if its
  \code{sparse} argument is false or true. In both cases, the resulting
  \code{"modelMatrix"} can then be coerced to a sparse or dense
  \code{"predModule"}.
}
\examples{
showClass("dPredModule")
showClass("sPredModule")

## see   example(model.Matrix)
}
\keyword{classes}
