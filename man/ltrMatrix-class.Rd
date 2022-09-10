\name{ltrMatrix-class}
\title{Triangular Dense Logical Matrices}
\docType{class}
%
\alias{ltpMatrix-class}
\alias{ltrMatrix-class}
%
\alias{!,ltpMatrix-method}
\alias{coerce,ltpMatrix,dtpMatrix-method}
\alias{coerce,ltpMatrix,lgeMatrix-method}
\alias{coerce,ltpMatrix,ltCMatrix-method}
\alias{coerce,ltpMatrix,ltRMatrix-method}
\alias{coerce,ltpMatrix,ltTMatrix-method}
\alias{coerce,ltpMatrix,lspMatrix-method}
\alias{coerce,ltpMatrix,ltrMatrix-method}
\alias{coerce,ltpMatrix,ntpMatrix-method}
\alias{coerce,matrix,ltpMatrix-method}
%
\alias{!,ltrMatrix-method}
\alias{coerce,ltrMatrix,dtrMatrix-method}
\alias{coerce,ltrMatrix,lgeMatrix-method}
\alias{coerce,ltrMatrix,ltCMatrix-method}
\alias{coerce,ltrMatrix,ltRMatrix-method}
\alias{coerce,ltrMatrix,ltTMatrix-method}
\alias{coerce,ltrMatrix,lsyMatrix-method}
\alias{coerce,ltrMatrix,ltpMatrix-method}
\alias{coerce,ltrMatrix,ntrMatrix-method}
\alias{coerce,matrix,ltrMatrix-method}
%
\description{
  The \code{"ltrMatrix"} class is the class of triangular, dense,
  logical matrices in nonpacked storage.  The \code{"ltpMatrix"} class
  is the same except in packed storage.
}
\section{Slots}{
  \describe{
    \item{\code{x}:}{Object of class \code{"logical"}. The logical
      values that constitute the matrix, stored in column-major order.}
    \item{\code{uplo}:}{Object of class \code{"character"}. Must be
      either "U", for upper triangular, and "L", for lower triangular.}
    \item{\code{diag}:}{Object of class \code{"character"}. Must be
      either \code{"U"}, for unit triangular (diagonal is all ones), or
      \code{"N"}; see \code{\linkS4class{triangularMatrix}}.}
    \item{\code{Dim},\code{Dimnames}:}{The dimension (a length-2
      \code{"integer"}) and corresponding names (or \code{NULL}), see the
      \code{\linkS4class{Matrix}} class.}
    \item{\code{factors}:}{Object of class \code{"list"}.  A named
      list of factorizations that have been computed for the matrix.}
  }
}
\section{Extends}{
  Both extend classes \code{"\linkS4class{ldenseMatrix}"} and
  \code{"\linkS4class{triangularMatrix}"}, directly; further, class
  \code{"Matrix"}, \code{"\linkS4class{lMatrix}"} and others,
  \emph{in}directly.  Use \code{\link{showClass}("ltrMatrix")}, e.g.,
  for details.
}
\section{Methods}{
  Currently, mainly \code{\link{t}()} and coercion methods (for
  \code{\link{as}(.)}; use, e.g.,
  \code{\link{showMethods}(class="ltpMatrix")} for details.
}
\seealso{
  Classes \code{\linkS4class{lgeMatrix}}, \code{\linkS4class{Matrix}};
  function \code{\link[base]{t}}
}
\examples{
showClass("ltrMatrix")

str(new("ltpMatrix"))
(lutr <- as(upper.tri(matrix(, 4, 4)), "ldenseMatrix"))
str(lutp <- pack(lutr)) # packed matrix: only 10 = 4*(4+1)/2 entries
!lutp # the logical negation (is *not* logical triangular !)
## but this one is:
stopifnot(all.equal(lutp, pack(!!lutp)))
}
\keyword{classes}
