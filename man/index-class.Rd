\name{index-class}
\title{Virtual Class \dQuote{index} of Index Vectors}
%
\docType{class}
\keyword{classes}
%
\alias{index-class}
%
\description{
  Class \code{index} is a virtual class designating index vectors,
  or \dQuote{subscripts}, for (possibly named) vectors and arrays.
  It is typically used in signatures of methods for the subscript
  and subassignment operators, namely \code{[} and \code{[<-}.
  It is implemented as a \emph{union} of the atomic vector classes
  \code{\linkS4class{numeric}}, \code{\linkS4class{logical}},
  and \code{\linkS4class{character}}.
}
\seealso{
  \code{\link{[}}, \code{\link{[-methods}}, and
  \code{\link{[<--methods}}.
}
\examples{
showClass("index")
}
