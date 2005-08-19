\name{Unused-classes}
%% This is just a stub, not yet used really: -- move to own page when used
\docType{class}
\alias{ldenseMatrix-class}
\alias{iMatrix-class}
\alias{zMatrix-class}
\title{Virtual Classes Not Yet Used}

\description{
  \code{ldenseMatrix} is the virtual class of all dense \bold{l}ogical
  (S4) matrices.  It extends both \code{\link{denseMatrix-class}}
  and \code{\link{lMatrix-class}} directly.

  \code{iMatrix} is the virtual class of all \bold{i}nteger
  (S4) matrices.  It extends the \code{\link{Matrix-class}} directly.

  \code{zMatrix} is the virtual class of all \code{\link{complex}}
  (S4) matrices.  It extends the \code{\link{Matrix-class}} directly.
}
\examples{
showClass("ldenseMatrix")

showClass("iMatrix")
showClass("zMatrix")
}
\keyword{classes}

