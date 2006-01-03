\name{ddiMatrix-class}
\docType{class}
\alias{ddiMatrix-class}
\alias{\%*\%,ddiMatrix,ddiMatrix-method}
\alias{crossprod,ddiMatrix,ddiMatrix-method}
\alias{tcrossprod,ddiMatrix,ddiMatrix-method}

\title{Class "ddiMatrix" of Diagonal Numeric Matrices}
\description{The class \code{"ddiMatrix"} of numerical diagonal matrices.
  %% FIXME add more
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("ddiMatrix", ...)}
  but typically rather via \code{\link{Diagonal}}.
}
\section{Slots}{
  \describe{
    \item{\code{x}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{diag}:}{Object of class \code{"character"} ~~ }
    \item{\code{Dim}:}{Object of class \code{"integer"} ~~ }
    \item{\code{Dimnames}:}{Object of class \code{"list"} ~~ }
    \item{\code{factors}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Extends}{
Class \code{"diagonalMatrix"}, directly.
Class \code{"dMatrix"}, directly.
Class \code{"denseMatrix"}, by class \code{"diagonalMatrix"}.
Class \code{"Matrix"}, by class \code{"denseMatrix"} and class \code{"dMatrix"}.
}
\section{Methods}{
  \describe{
    \item{\%*\%}{\code{signature(x = "ddiMatrix", y = "ddiMatrix")}: ... }
	 }
}
\seealso{
  \code{\link{diagonalMatrix-class}};
  \code{\link{Diagonal}}
}
\examples{
(d2 <- Diagonal(x = c(10,1)))
str(d2)

M <- Matrix(cbind(1,2:4))
M \%*\% d2 #> `fast' multiplication
}
\keyword{classes}
