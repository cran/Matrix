\name{ldiMatrix-class}
\docType{class}
\alias{ldiMatrix-class}
\alias{!,ldiMatrix-method}

\title{Class "ldiMatrix" of Diagonal Logical Matrices}
\description{The class \code{"ldiMatrix"} of logical diagonal matrices.
%% FIXME add more
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("ldiMatrix", ...)}
  but typically rather via \code{\link{Diagonal}}.
}
\section{Slots}{
	 \describe{
    \item{\code{x}:}{Object of class \code{"logical"} ~~ }
    \item{\code{diag}:}{Object of class \code{"character"} ~~ }
    \item{\code{Dim}:}{Object of class \code{"integer"} ~~ }
    \item{\code{Dimnames}:}{Object of class \code{"list"} ~~ }
    \item{\code{factors}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Extends}{
Class \code{"diagonalMatrix"}, directly.
Class \code{"lMatrix"}, directly.
Class \code{"denseMatrix"}, by class \code{"diagonalMatrix"}.
Class \code{"Matrix"}, by class \code{"denseMatrix"} and class \code{"lMatrix"}.
}
% \section{Methods}{
% No methods defined with class "ldiMatrix" in the signature.
% }
% \references{ ~put references to the literature/web site here ~ }
\seealso{
  \code{\link{ddiMatrix-class}},
  \code{\link{diagonalMatrix-class}};
  \code{\link{Diagonal}}
}
\examples{
(lM <- Diagonal(x = c(TRUE,FALSE,FALSE)))
str(lM)#> gory details (slots)
%% FIXME
#not yet: crossprod(lM)
}
\keyword{classes}
