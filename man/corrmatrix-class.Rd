\name{corrmatrix-class}
\docType{class}
\alias{corrmatrix-class}
\title{Class "corrmatrix"}
\alias{coerce,pdIdent,corrmatrix-method}
\alias{coerce,pdMat,corrmatrix-method}
\alias{coerce,pdNatural,corrmatrix-method}
\alias{coerce<-,pdNatural,corrmatrix-method}
\alias{coerce,pdmatrix,corrmatrix-method}
\description{A class of correlation matrices}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("corrmatrix", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"matrix"} - the actual
      correlation matrix}
    \item{\code{stdDev}:}{Object of class \code{"numeric"} giving the
      corresponding standard deviations }
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
    \item{coerce}{\code{signature(from = "pdIdent", to = "corrmatrix")}: ... }
  }
}
\references{}
\author{Douglas Bates \email{bates@stat.wisc.edu}
  and Saikat DebRoy \email{saikat@stat.wisc.edu}}
\examples{
}
\keyword{classes}
