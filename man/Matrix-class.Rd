\name{Matrix-class}
\docType{class}
\alias{Matrix-class}
\alias{-,Matrix,missing-method}
\alias{+,Matrix,missing-method}
\alias{\%*\%,Matrix,numeric-method}
\alias{\%*\%,numeric,Matrix-method}
\alias{\%*\%,Matrix,Matrix-method}
% "[" are in ./Xtrct-methods.Rd <<<<<<<<
\alias{cbind2,Matrix,Matrix-method}
\alias{cbind2,Matrix,NULL-method}
\alias{cbind2,Matrix,missing-method}
\alias{cbind2,NULL,Matrix-method}
\alias{rbind2,Matrix,Matrix-method}
\alias{rbind2,Matrix,NULL-method}
\alias{rbind2,Matrix,missing-method}
\alias{rbind2,NULL,Matrix-method}
\alias{crossprod,Matrix,numeric-method}
\alias{crossprod,Matrix,ANY-method}
\alias{crossprod,numeric,Matrix-method}
\alias{solve,Matrix,numeric-method}
\alias{dim,Matrix-method}
\alias{dimnames,Matrix-method}
\alias{dimnames<-,Matrix,list-method}
\alias{show,Matrix-method}
\alias{as.array,Matrix-method}
\alias{as.matrix,Matrix-method}
\alias{as.vector,Matrix,missing-method}
\alias{t,Matrix-method}
\alias{unname,Matrix-method}
%\alias{coerce,Matrix,matrix-method}
\alias{Compare,Matrix,Matrix-method}
\alias{Compare,Matrix,ANY-method}
\alias{Compare,ANY,Matrix-method}
\alias{crossprod,ANY,Matrix-method}
\alias{kronecker,ANY,Matrix-method}
\alias{kronecker,Matrix,ANY-method}
\alias{[,Matrix,lMatrix,missing,ANY-method}
\alias{[,Matrix,logical,missing,ANY-method}
%
\title{Virtual Class "Matrix" Class of Matrices}
\description{
  The \code{Matrix} class is a class contained by all actual
  classes in the \pkg{Matrix} package.  It is a \dQuote{virtual} class.
}
\section{Slots}{
  Common to \emph{all} matrix objects in the package:
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
    \item{dim}{\code{signature(x = "Matrix")}: extract matrix dimensions
      \code{\link{dim}}.}
    \item{dimnames}{\code{signature(x = "Matrix")}: extract
      \code{\link{dimnames}}.}
    \item{dimnames<-}{\code{signature(x = "Matrix", value = "list")}: set
      the \code{dimnames} to a \code{\link{list}} of length 2, see
      \code{\link{dimnames<-}}.}
    \item{show}{\code{signature(object = "Matrix")}: \code{\link{show}}
      method for \code{\link{print}}ing.}
    \cr %------------------------------------
    \item{as.matrix, as.array}{\code{signature(x = "Matrix")}: the same as
      \code{as(x, "matrix")}; see also the note below.}
    \item{as.vector}{\code{signature(x = "Matrix", mode = "missing")}:
      \code{as.vector(m)} should be identical to \code{as.vector(as(m,
	"matrix"))}, implemented more efficiently for some subclasses.}
  }
}
\note{
  Loading the \code{Matrix} namespace \dQuote{overloads}
  \code{\link{as.matrix}} and \code{\link{as.array}} in the \pkg{base}
  namespace by the equivalent of \code{function(x) as(x, "matrix")}.
  Consequently, \code{as.matrix(m)} or \code{as.array(m)} will properly
  work when \code{m} inherits from the \code{"Matrix"} class ---
  \emph{also} for functions in package \pkg{base} and other packages.
  E.g., \code{\link{apply}} or \code{\link{outer}} can therefore be applied
  to \code{"Matrix"} matrices.
}
%\references{}
\author{Douglas Bates \email{bates@stat.wisc.edu} and Martin Maechler}
\seealso{
  \code{\link{dgeMatrix-class}}, \code{\link{dgCMatrix-class}}, and
  \code{\link{Matrix}} for construction (and examples).
}
\examples{
slotNames("Matrix")

cl <- getClass("Matrix")
names(cl@subclasses) # more than 40 ..

showClass("Matrix")#> output with slots and all subclasses
}
\keyword{classes}
\keyword{algebra}
