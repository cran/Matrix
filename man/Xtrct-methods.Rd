\name{[-methods}
\docType{methods}
\alias{[-methods}
\title{Methods for "[": Extraction or Subsetting in Package 'Matrix'}
%
\alias{[,Matrix,missing,missing,ANY-method}
\alias{[,Matrix,missing,numeric,missing-method}
\alias{[,Matrix,numeric,missing,missing-method}
\alias{[,Matrix,numeric,numeric,missing-method}
%
\alias{[,gTMatrix,numeric,missing,logical-method}
\alias{[,gTMatrix,missing,numeric,logical-method}
\alias{[,gTMatrix,numeric,numeric,logical-method}
%
\alias{[,dgTMatrix,missing,missing,ANY-method}
\alias{[,dgTMatrix,missing,numeric,logical-method}
\alias{[,dgTMatrix,missing,numeric,missing-method}
\alias{[,dgTMatrix,numeric,missing,logical-method}
\alias{[,dgTMatrix,numeric,missing,missing-method}
\alias{[,dgTMatrix,numeric,numeric,logical-method}
\alias{[,dgTMatrix,numeric,numeric,missing-method}
%
\alias{[,dsparseMatrix,missing,numeric,logical-method}
\alias{[,dsparseMatrix,numeric,missing,logical-method}
\alias{[,dsparseMatrix,numeric,numeric,logical-method}
%
\alias{[,lsparseMatrix,missing,numeric,logical-method}
\alias{[,lsparseMatrix,numeric,missing,logical-method}
\alias{[,lsparseMatrix,numeric,numeric,logical-method}
%-------
\description{
  Methods for \code{"["}, i.e., extraction or subsetting mostly of
  matrices, in package \pkg{Matrix}.
}
\section{Methods}{
  There are more than these:
  \describe{
    \item{x = "Matrix", i = "missing", j = "missing", drop= "ANY"}{ ... }
    \item{x = "Matrix", i = "numeric", j = "missing", drop= "missing"}{ ... }
    \item{x = "Matrix", i = "missing", j = "numeric", drop= "missing"}{ ... }

    \item{x = "dsparseMatrix", i = "missing", j = "numeric", drop= "logical"}{ ... }
    \item{x = "dsparseMatrix", i = "numeric", j = "missing", drop= "logical"}{ ... }
    \item{x = "dsparseMatrix", i = "numeric", j = "numeric", drop= "logical"}{ ... }
  }
}
\seealso{
  \code{\link{Extract}} about the standard extraction.
}
\examples{
str(m <- Matrix(round(rnorm(7*4),2), nrow = 7))
stopifnot(identical(m, m[]))
\dontshow{
 if(FALSE) {## FIXME: still all fail
   m[2,3]
   m[1,]
   m[,2]
 }
}
}
\keyword{methods}
\keyword{array}
