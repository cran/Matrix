\name{isInitialized}
\docType{genericFunction}
\alias{isInitialized}
\alias{isInitialized-methods}
\alias{isInitialized,pdBlocked-method}
\alias{isInitialized,pdCompSymm-method}
\alias{isInitialized,pdMat-method}
\alias{isInitialized,pdNatural-method}
%\alias{isInitialized,pdScalar-method}
\alias{isInitialized,lmeLevelList-method}
\title{Check for initialization}
\description{
  A generic function.  Methods for specific classes should return
  \code{TRUE} if \code{object} is initialized, otherwise \code{FALSE}.
}
\usage{
isInitialized(object)
}
\arguments{
  \item{object}{An object that may or may not be initialized.}
}
\value{
  A logical value indicating if the object is initialized.
}
\note{The name of this generic was changed from the original name,
  \code{isinitialized}, because that name is used in the \code{methods}
  package. }
\references{}
\author{Douglas Bates \email{bates@stat.wisc.edu}
  and Saikat DebRoy \email{saikat@stat.wisc.edu}}
\examples{
}
\keyword{manip}
