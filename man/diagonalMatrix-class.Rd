\name{diagonalMatrix-class}
\docType{class}
\alias{diagonalMatrix-class}
\alias{coerce,matrix,diagonalMatrix-method}
\alias{coerce,diagonalMatrix,generalMatrix-method}
\alias{coerce,diagonalMatrix,triangularMatrix-method}
\alias{coerce,diagonalMatrix,matrix-method}
\alias{coerce,Matrix,diagonalMatrix-method}
\alias{t,diagonalMatrix-method}
\alias{diag,diagonalMatrix-method}
\alias{show,diagonalMatrix-method}
\alias{\%*\%,dgeMatrix,diagonalMatrix-method}
\alias{\%*\%,matrix,diagonalMatrix-method}
\alias{\%*\%,diagonalMatrix,matrix-method}
\alias{\%*\%,diagonalMatrix,dgeMatrix-method}
\alias{crossprod,diagonalMatrix,matrix-method}
\alias{crossprod,diagonalMatrix,dgeMatrix-method}

\title{Class "diagonalMatrix" of Diagonal Matrices}
\description{
  Class "diagonalMatrix" is the virtual class of all diagonal matrices.
}
\section{Objects from the Class}{A virtual Class: No objects may be
  created from it.}
\section{Slots}{
  \describe{
    \item{\code{diag}:}{Object of class \code{"character"} ~~ }
    \item{\code{Dim}:}{Object of class \code{"integer"} ~~ }
    \item{\code{Dimnames}:}{a \code{\link{list}}, typically
      \code{list(NULL,NULL)} for diagonal matrices.}
    \item{\code{factors}:}{a \code{\link{list}}, typically empty.}
  }
}
\section{Extends}{
  Class \code{"denseMatrix"}, directly.
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "matrix", to = "diagonalMatrix")}: ... }
    \item{coerce}{\code{signature(from = "diagonalMatrix", to = "triangularMatrix")}: ... }
    \item{coerce}{\code{signature(from = "diagonalMatrix", to = "matrix")}: ... }
    \item{coerce}{\code{signature(from = "Matrix", to = "diagonalMatrix")}: ... }
    \item{t}{\code{signature(x = "diagonalMatrix")}: ... }
    \item{\%*\%}{\code{signature(x = "dgeMatrix", y = "diagonalMatrix")}: ... }
    \item{\%*\%}{\code{signature(x = "matrix", y = "diagonalMatrix")}: ... }
    \item{\%*\%}{\code{signature(x = "diagonalMatrix", y = "matrix")}: ... }
    \item{\%*\%}{\code{signature(x = "diagonalMatrix", y = "dgeMatrix")}: ... }
	 }
}
\seealso{
  \code{\link{ddiMatrix-class}} and \code{\link{ldiMatrix-class}} are
  \dQuote{actual} classes extending \code{"diagonalMatrix"}
}
% \examples{
% }
\keyword{classes}
