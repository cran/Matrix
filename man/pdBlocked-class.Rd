\name{pdBlocked-class}
\docType{class}
\alias{pdBlocked-class}
\alias{pdBlocked}
\title{Class "pdBlocked", blocked pd matrices}
\description{
  Block-diagonal positive definite matrices.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("pdBlocked", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{param}:}{Object of class \code{"list"} }
    \item{\code{form}:}{Object of class \code{"formula", from class "pdMat"}} \item{\code{Names}:}{Object of class \code{"character", from class "pdMat"} }
  }
}
\section{Extends}{
  Class \code{"pdMat"}, directly.
}
\section{Methods}{
  \describe{
    \item{Names}{\code{signature(object = "pdBlocked")}:}
    \item{Names<-}{\code{signature(object = "pdBlocked", value = "character")}:}
    \item{[}{\code{signature(x = "pdBlocked")}: }
    \item{coef}{\code{signature(object = "pdBlocked")}: }
    \item{coef<-}{\code{signature(object = "pdBlocked", value = "ANY")}: }
    \item{corMatrix}{\code{signature(object = "pdBlocked")}:}
    \item{formula}{\code{signature(x = "pdBlocked")}:}
    \item{isInitialized}{\code{signature(object = "pdBlocked")}:}
    \item{logDet}{\code{signature(object = "pdBlocked", covariate = "missing")}:}
    \item{matrix<-}{\code{signature(object = "pdBlocked", value = "matrix")}: }
    \item{pdFactor}{\code{signature(object = "pdBlocked")}: }
    \item{solve}{\code{signature(a = "pdBlocked", b = "missing")}:}
    \item{summary}{\code{signature(object = "pdBlocked")}:}
  }
}
%\references{}
\author{Douglas Bates \email{bates@stat.wisc.edu}
  and Saikat DebRoy \email{saikat@stat.wisc.edu}}
%\seealso{}
%\examples{}
\keyword{classes}
