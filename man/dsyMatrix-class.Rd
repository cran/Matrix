\name{dsyMatrix-class}
\docType{class}
\alias{\%*\%,dgeMatrix,dsyMatrix-method}
\alias{\%*\%,dspMatrix,dgeMatrix-method}
\alias{\%*\%,dspMatrix,integer-method}
\alias{\%*\%,dspMatrix,matrix-method}
\alias{\%*\%,dspMatrix,numeric-method}
\alias{\%*\%,dsyMatrix,dgeMatrix-method}
\alias{coerce,dspMatrix,dgeMatrix-method}
\alias{coerce,dspMatrix,dsyMatrix-method}
\alias{coerce,dspMatrix,lspMatrix-method}
\alias{coerce,dspMatrix,matrix-method}
\alias{coerce,dgeMatrix,dspMatrix-method}
\alias{coerce,matrix,dspMatrix-method}
\alias{coerce,dgeMatrix,dsyMatrix-method}% or rather setIs?
\alias{coerce,dsyMatrix,dgeMatrix-method}
\alias{coerce,dsyMatrix,dspMatrix-method}
\alias{coerce,dsyMatrix,lsyMatrix-method}
\alias{coerce,dsyMatrix,matrix-method}
\alias{coerce,matrix,dsyMatrix-method}
\alias{show,dsyMatrix-method}
\alias{dspMatrix-class}
\alias{dsyMatrix-class}
\alias{norm,dspMatrix,character-method}
\alias{norm,dspMatrix,missing-method}
\alias{norm,dsyMatrix,character-method}
\alias{norm,dsyMatrix,missing-method}
\alias{rcond,dspMatrix,character-method}
\alias{rcond,dspMatrix,missing-method}
\alias{rcond,dsyMatrix,character-method}
\alias{rcond,dsyMatrix,missing-method}
\alias{solve,dspMatrix,dgeMatrix-method}
\alias{solve,dspMatrix,integer-method}
\alias{solve,dspMatrix,matrix-method}
\alias{solve,dspMatrix,missing-method}
\alias{solve,dspMatrix,numeric-method}
\alias{solve,dsyMatrix,dgeMatrix-method}
\alias{solve,dsyMatrix,matrix-method}
\alias{solve,dsyMatrix,missing-method}
\alias{t,dspMatrix-method}
\alias{t,dsyMatrix-method}
\alias{unpack,dspMatrix-method}
\title{Symmetric Dense Numeric Matrices}
\description{
  The \code{"dsyMatrix"} class is the class of symmetric, dense matrices
  in non-packed storage and \code{"dspMatrix"} is the class of symmetric
  dense matrices in packed storage.  Only the upper triangle or the
  lower triangle is stored.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("dsyMatrix", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{uplo}:}{Object of class \code{"character"}. Must be
      either "U", for upper triangular, and "L", for lower triangular.}
    \item{\code{x}:}{Object of class \code{"numeric"}. The numeric
      values that constitute the matrix, stored in column-major order.}
    \item{\code{Dim},\code{Dimnames}:}{The dimension (a length-2
      \code{"integer"}) and corresponding names (or \code{NULL}), see the
      \code{\linkS4class{Matrix}}.}
    \item{\code{factors}:}{Object of class \code{"list"}.  A named
      list of factorizations that have been computed for the matrix.}
  }
}
\section{Extends}{
  \code{"dsyMatrix"} extends class \code{"dgeMatrix"}, directly, whereas\cr
  \code{"dspMatrix"} extends class \code{"ddenseMatrix"}, directly.

  Both extend class \code{"symmetricMatrix"}, directly,
  and class \code{"Matrix"} and others, \emph{in}directly, use
  \code{\link{showClass}("dsyMatrix")}, e.g., for details.
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "dspMatrix", to = "dgeMatrix")}}
    \item{coerce}{\code{signature(from = "dsyMatrix", to = "dgeMatrix")}}
    \item{coerce}{\code{signature(from = "dspMatrix", to = "matrix")}}
    \item{coerce}{\code{signature(from = "dsyMatrix", to = "matrix")}}
    \item{coerce}{\code{signature(from = "dsyMatrix", to = "dspMatrix")}}
    \item{coerce}{\code{signature(from = "dspMatrix", to = "dsyMatrix")}}
    \item{norm}{\code{signature(x = "dspMatrix", type = "character")}}
    \item{norm}{\code{signature(x = "dsyMatrix", type = "character")}}
    \item{norm}{\code{signature(x = "dspMatrix", type = "missing")}}
    \item{norm}{\code{signature(x = "dsyMatrix", type = "missing")}}
    \item{solve}{\code{signature(a = "dspMatrix", b = "missing")}}
    \item{solve}{\code{signature(a = "dsyMatrix", b = "missing")}}
    \item{solve}{\code{signature(a = "dspMatrix", b = "dgeMatrix")}}
    \item{solve}{\code{signature(a = "dsyMatrix", b = "dgeMatrix")}}
    \item{solve}{\code{signature(a = "dspMatrix", b = "matrix")}}
    \item{solve}{\code{signature(a = "dsyMatrix", b = "matrix")}}
  }
}
%\references{}
%\author{}
\seealso{
  Classes \code{\linkS4class{dgeMatrix}} and \code{\linkS4class{Matrix}};
  \code{\link[base]{solve}}, \code{\link{norm}}, \code{\link{rcond}},
  \code{\link[base]{t}}
}
\examples{
## Only upper triangular part matters (when uplo == "U" as per default)
(sy2 <- new("dsyMatrix", Dim = as.integer(c(2,2)), x = c(14, NA,32,77)))
(c2 <- chol(sy2))
all.equal(c2, c2. <- chol(as(sy2, "dpoMatrix")))# equal, but not same class
str(c2)

## An example where chol() can't work
(sy3 <- new("dsyMatrix", Dim = as.integer(c(2,2)), x = c(14, -1, 2, -7)))

validObject(as(sy3, "dpoMatrix"), test=TRUE) # >> is not pos.def.
% .... Lapack routine dpotrf returned error code 2 -- FIXME: better message
try(chol(sy3)) ## Error: not pos.def
}
\keyword{classes}
