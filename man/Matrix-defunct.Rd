\name{Matrix-defunct}
\title{Defunct Functions in Package \pkg{Matrix}}
%
\alias{Matrix-defunct}
%
\alias{cBind}
\alias{rBind}
%
\description{
  The functions or variables listed here are no longer part of
  \pkg{Matrix} as they are no longer needed.
}
\usage{
## Defunct in 1.3-3
rBind(..., deparse.level = 1)
cBind(..., deparse.level = 1)
}
\details{
  These either are stubs reporting that they are defunct, or
  have been removed completely (apart from being documented here).

  \code{rBind} and \code{cBind} were provided for \R{} versions
  older than 3.2.0 as substitutes for the \pkg{base} analogues
  \code{rbind} and \code{cbind}, which at that time were not
  \dQuote{S4-aware} and so could not be used to vertically or
  horizontally concatenate \code{"\linkS4class{Matrix}"} objects
  together with traditional matrices and vectors.
}
\seealso{
  \code{\link{Defunct}},
  \code{\link{base-defunct}},
  \code{\link{Matrix-deprecated}}
}
\keyword{internal}
