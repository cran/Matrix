\name{modelMatrix-class}
\Rdversion{1.1}
\title{Class "modelMatrix" and SubClasses}
\docType{class}
\alias{modelMatrix-class}
\alias{denseModelMatrix-class}
\alias{ddenseModelMatrix-class}
\alias{sparseModelMatrix-class}
\alias{dsparseModelMatrix-class}
%
\alias{show,modelMatrix-method}
\description{
  The class \code{"modelMatrix"} and notably its subclass
  \code{"dsparseModelMatrix"} are used to encode additional information,
  analogously to what the standard \R function
  \code{\link{model.matrix}()} returns.
}
\section{Objects from the Classes}{
  Only \code{"dsparseModelMatrix"} and \code{"ddenseModelMatrix"} are
  \dQuote{actual} (aka non-virtual) classes. For these,
  objects can be created by calls of the form
  \code{new("dsparseModelMatrix", x, assign, contrast)}, where \code{x}
  is a \code{\linkS4class{dgCMatrix}} classed object.
}
\section{Slots}{
  The \code{"modelMatrix"} mother class contains
  \code{\linkS4class{Matrix}} plus two extra slots,
  \describe{
    \item{\code{assign}:}{\code{"integer"} vector of length
      \code{ncol(.)}, coding the variables which make up the
      matrix columns, see \code{\link{model.matrix}}.}
    \item{\code{contrasts}:}{a named \code{\link{list}} of
      \code{\link{contrasts}}, as in \code{\link{model.matrix}()}.}

    \item{\code{Dim}:}{integer vector of length two with the matrix dimensions.}
    \item{\code{Dimnames}:}{list of length two, the
      \code{\link{dimnames}(.)} of the matrix.}
  }
  whereas the (current only) actual classes \code{"d*ModelMatrix"},
  have an at least an additional (\code{\link{numeric}} slot \code{"x"}.
  E.g., "dsparseModelMatrix" has the additional slots
  \describe{
    \item{\code{i},\code{p}:}{row number and \dQuote{pointer} integer
      vectors, see class \code{"\linkS4class{dgCMatrix}"}.}
    \item{\code{x}:}{\code{"numeric"} vector of non-zero entries.}
    \item{\code{factors}:}{a (possibly empty) \code{\link{list}} of
      factorizations.}
  }
}
\section{Extends}{
  \code{"dsparseModelMatrix"} extends class \code{"\linkS4class{dgCMatrix}"} directly,\cr
   \code{"ddenseModelMatrix"} extends class \code{"\linkS4class{dgeMatrix}"} directly.
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "modelMatrix")}:
      \code{\link{show}(.)} the matrix, but also the \code{assign} and
      \code{contrasts} slots.}
  }
}
\author{Martin Maechler}
\seealso{
  \code{\link{sparse.model.matrix}} will return a
  \code{"dsparseModelMatrix"} object.
  \code{\link{model.Matrix}} which is a simple wrapper around the traditional
  \code{\link{model.matrix}} and returns a \code{"ddenseModelMatrix"} object.
}
\examples{
showClass("modelMatrix")
showClass("dsparseModelMatrix")

## see   example(model.Matrix)
}
\keyword{classes}
