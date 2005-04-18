\name{[-methods}
\docType{methods}
\alias{[-methods}
%
\alias{[,dgTMatrix,missing,missing,ANY-method}
\alias{[,dgTMatrix,missing,numeric,logical-method}
\alias{[,dgTMatrix,missing,numeric,missing-method}
\alias{[,dgTMatrix,numeric,missing,logical-method}
\alias{[,dgTMatrix,numeric,missing,missing-method}
\alias{[,dgTMatrix,numeric,numeric,logical-method}
\alias{[,dgTMatrix,numeric,numeric,missing-method}
%
\alias{[,dsparseMatrix,missing,missing,ANY-method}
\alias{[,dsparseMatrix,missing,numeric,logical-method}
\alias{[,dsparseMatrix,missing,numeric,missing-method}
\alias{[,dsparseMatrix,numeric,missing,logical-method}
\alias{[,dsparseMatrix,numeric,missing,missing-method}
\alias{[,dsparseMatrix,numeric,numeric,logical-method}
\alias{[,dsparseMatrix,numeric,numeric,missing-method}
\title{Methods for "[": Extraction or Subsetting in Package 'Matrix'.
}
\description{
  Methods for \code{"["}, i.e., extraction or subsetting in package
  \pkg{Matrix}.
}
\section{Methods}{
  There are more than these:
  \describe{
    \item{x = "dsparseMatrix", i = "missing", j = "missing", drop= "ANY"}{ ... }
    \item{x = "dsparseMatrix", i = "missing", j = "numeric", drop= "logical"}{ ... }
    \item{x = "dsparseMatrix", i = "missing", j = "numeric", drop= "missing"}{ ... }
    \item{x = "dsparseMatrix", i = "numeric", j = "missing", drop= "logical"}{ ... }
    \item{x = "dsparseMatrix", i = "numeric", j = "missing", drop= "missing"}{ ... }

    \item{x = "dsparseMatrix", i = "numeric", j = "numeric", drop= "logical"}{ ... }
}}
\seealso{
  \code{\link{Extract}} about the standard extraction.
}
\examples{
str(m <- Matrix(round(rnorm(7*4),2), nrow = 7))
## wrongly looses class: m[]; error: m[1,1]

\dontshow{
 if(FALSE) {## still all fail
   ch[2,3]
   ch[1,]
   ch[,2]
   pm[2,3]
   pm[1,]
   pm[,2]
 }
}
}
\keyword{methods}
\keyword{array}
