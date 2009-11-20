\name{abIndex-class}
\title{Class "abIndex" of Abstract Index Vectors}
\Rdversion{1.1}
\docType{class}
%
\alias{abIndex-class}
%
\alias{as.numeric,abIndex-method}
\alias{[,abIndex,index,ANY,ANY-method}
\alias{coerce,numeric,abIndex-method}
\alias{coerce,abIndex,numeric-method}
\alias{coerce,abIndex,integer-method}
\alias{length,abIndex-method}
\alias{Ops,numeric,abIndex-method}
\alias{Ops,abIndex,abIndex-method}
\alias{Ops,abIndex,numeric-method}
\alias{Summary,abIndex-method}
%
\description{
  The \code{"abIndex"} \code{\link{class}}, short for \dQuote{Abstract
    Index Vector}, is used for dealing with large index vectors more
  efficiently, than using integer (or \code{\link{numeric}}) vectors of
  the kind \code{2:1000000} or \code{c(0:1e5, 1000:1e6)}.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("abIndex", ...)},
  or more easily by \code{as(x, "abIndex")} where \code{x} is an integer
  (valued) vector.
}
\section{Slots}{
  \describe{
    \item{\code{kind}:}{Object of class \code{"character"} ~~ }
    \item{\code{x}:}{Object of class \code{"numIndex"} ~~ }
    \item{\code{rleD}:}{Object of class \code{"rleDiff"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{as.numeric}{\code{signature(x = "abIndex")}: ... }
    \item{[}{\code{signature(x = "abIndex", i = "index", j = "ANY", drop = "ANY")}: ... }
    \item{coerce}{\code{signature(from = "numeric", to = "abIndex")}: ... }
    \item{coerce}{\code{signature(from = "abIndex", to = "numeric")}: ... }
    \item{coerce}{\code{signature(from = "abIndex", to = "integer")}: ... }
    \item{length}{\code{signature(x = "abIndex")}: ... }
    \item{Ops}{\code{signature(e1 = "numeric", e2 = "abIndex")}: These
      and the following arithmetic and logic operations are
      \bold{not yet implemented}; see \code{\link[methods]{Ops}} for a
      list of these (S4) group methods.}
    \item{Ops}{\code{signature(e1 = "abIndex", e2 = "abIndex")}: ... }
    \item{Ops}{\code{signature(e1 = "abIndex", e2 = "numeric")}: ... }
    \item{Summary}{\code{signature(x = "abIndex")}: ... }
	 }
}
\note{
  This is currently experimental and not yet used for our own code.
  Please contact us (\code{packageDescription("Matrix")$Maintainer}),
  if you plan to make use of this class.

  Partly builds on ideas and code from Jens Oehlschlaegel,
  as implemented (around 2008, in the GPL'ed part of) package \pkg{ff}.
}
%\author{Martin Maechler}
\seealso{
  \code{\link{rle}} (\pkg{base}) which is used here.
  \code{\link{numeric}}
	% ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
	% or \code{\linkS4class{CLASSNAME}} for links to other classes
}
\examples{
showClass("abIndex")
ii <- c(-3:40, 20:70)
str(ai <- as(ii, "abIndex"))# note

%% FIXME: add / exchange with ../tests/abIndex-tsts.R
}
\keyword{classes}
