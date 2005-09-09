\name{Matrix-class}
\docType{class}
\alias{Matrix-class}
\alias{\%*\%,Matrix,numeric-method}
\alias{\%*\%,numeric,Matrix-method}
% "[" are in ./Xtrct-methods.Rd
\alias{cbind2,Matrix,Matrix-method}
\alias{cbind2,Matrix,NULL-method}
\alias{cbind2,Matrix,missing-method}
\alias{cbind2,NULL,Matrix-method}
\alias{crossprod,Matrix,numeric-method}
\alias{crossprod,numeric,Matrix-method}
\alias{solve,Matrix,numeric-method}
\alias{dim,Matrix-method}
\alias{dimnames,Matrix-method}
\alias{dimnames<-,Matrix,list-method}
\alias{show,Matrix-method}
\alias{as.matrix,Matrix-method}
\alias{unname,Matrix-method}
%\alias{coerce,Matrix,matrix-method}
\title{Virtual Class "Matrix" Class of Matrices}
\description{
  The \code{Matrix} class is a class contained by all actual
  classes in the \pkg{Matrix} package.  It is a \dQuote{virtual} class.
}
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
    \item{dim}{\code{signature(x = "Matrix")}: extract matrix dimensions
      \code{\link{dim}}.}
    \item{dimnames}{\code{signature(x = "Matrix")}: extract
      \code{\link{dimnames}}.}
    \item{dimnames<-}{\code{signature(x = "Matrix", value = "list")}: set
      the \code{dimnames} to a \code{\link{list}} of length 2, see
      \code{\link{dimnames<-}}.}
    \item{show}{\code{signature(object = "Matrix")}: \code{\link{show}}
      method for \code{\link{print}}ing.}
  }
}
%\references{}
\author{Douglas Bates \email{bates@stat.wisc.edu}}
\seealso{
  \code{\link{dgeMatrix-class}}, \code{\link{dgCMatrix-class}}, and
  \code{\link{Matrix}} for construction (and examples).
}
\keyword{classes}
\keyword{algebra}
