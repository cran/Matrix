\name{ldiMatrix-class}
\title{Class "ldiMatrix" of Diagonal Logical Matrices}
%
\docType{class}
\keyword{array}
\keyword{classes}
%
\alias{ldiMatrix-class}
%
\alias{!,ldiMatrix-method}
\alias{\%\%,ldiMatrix,Matrix-method}
\alias{\%\%,ldiMatrix,ddenseMatrix-method}
\alias{\%\%,ldiMatrix,ldenseMatrix-method}
\alias{\%\%,ldiMatrix,ndenseMatrix-method}
\alias{\%/\%,ldiMatrix,Matrix-method}
\alias{\%/\%,ldiMatrix,ddenseMatrix-method}
\alias{\%/\%,ldiMatrix,ldenseMatrix-method}
\alias{\%/\%,ldiMatrix,ndenseMatrix-method}
\alias{&,ldiMatrix,Matrix-method}
\alias{&,ldiMatrix,ddenseMatrix-method}
\alias{&,ldiMatrix,ldenseMatrix-method}
\alias{&,ldiMatrix,ndenseMatrix-method}
\alias{*,ldiMatrix,Matrix-method}
\alias{*,ldiMatrix,ddenseMatrix-method}
\alias{*,ldiMatrix,ldenseMatrix-method}
\alias{*,ldiMatrix,ndenseMatrix-method}
\alias{-,ldiMatrix,missing-method}
\alias{/,ldiMatrix,Matrix-method}
\alias{/,ldiMatrix,ddenseMatrix-method}
\alias{/,ldiMatrix,ldenseMatrix-method}
\alias{/,ldiMatrix,ndenseMatrix-method}
\alias{Arith,ldiMatrix,logical-method}
\alias{Arith,ldiMatrix,numeric-method}
\alias{Arith,logical,ldiMatrix-method}
\alias{Arith,numeric,ldiMatrix-method}
\alias{Ops,ANY,ldiMatrix-method}
\alias{Ops,ldiMatrix,ANY-method}
\alias{Ops,ldiMatrix,Matrix-method}
\alias{Ops,ldiMatrix,dMatrix-method}
\alias{Ops,ldiMatrix,ddiMatrix-method}
\alias{Ops,ldiMatrix,ldiMatrix-method}
\alias{Ops,ldiMatrix,logical-method}
\alias{Ops,ldiMatrix,numeric-method}
\alias{Ops,ldiMatrix,sparseMatrix-method}
\alias{Summary,ldiMatrix-method}
\alias{prod,ldiMatrix-method}
\alias{sum,ldiMatrix-method}
\alias{which,ldiMatrix-method}
%
\description{The class \code{"ldiMatrix"} of logical diagonal matrices.
%% FIXME add more
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("ldiMatrix", ...)}
  but typically rather via \code{\link{Diagonal}}.
}
\section{Slots}{
  \describe{
    \item{\code{x}:}{\code{"logical"} vector.}
    \item{\code{diag}:}{\code{"character"} string, either "U" or "N",
      see \code{\linkS4class{ddiMatrix}}.}
    \item{\code{Dim},\code{Dimnames}:}{matrix dimension and
      \code{\link{dimnames}}, see the \code{\linkS4class{Matrix}} class
      description.}
  }
}
\section{Extends}{
  Class \code{"\linkS4class{diagonalMatrix}"} and
  class \code{"\linkS4class{lMatrix}"}, directly.
  
  Class \code{"\linkS4class{sparseMatrix}"}, by class \code{"diagonalMatrix"}.
}
% \section{Methods}{
% No methods defined with class "ldiMatrix" in the signature.
% }
% \references{ ~put references to the literature/web site here ~ }
\seealso{
  Classes \code{\linkS4class{ddiMatrix}} and
  \code{\linkS4class{diagonalMatrix}}; function \code{\link{Diagonal}}.
}
\examples{
\dontshow{ % for R_DEFAULT_PACKAGES=NULL
library(utils, pos = "package:base", verbose = FALSE)
}
(lM <- Diagonal(x = c(TRUE,FALSE,FALSE)))
str(lM)#> gory details (slots)

crossprod(lM) # numeric
(nM <- as(lM, "nMatrix"))# -> sparse (not formally ``diagonal'')
crossprod(nM) # logical sparse
}
