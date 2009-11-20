\name{rleDiff-class}
\Rdversion{1.1}
\docType{class}
\alias{rleDiff-class}
\title{Class "rleDiff" of  rle(diff(.))  Stored Vectors}
\description{
  Class \code{"rleDiff"} is for compactly storing long vectors
  which mainly consist of \emph{linear} stretches.  For such a vector
  \code{x}, \code{\link{diff}(x)} consists of \emph{constant} stretches
  and is hence well compressable via \code{\link{rle}()}.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("rleDiff", ...)}.

  Currently experimental, see below.
}
\section{Slots}{
  \describe{
    \item{\code{first}:}{A single number (of class \code{"numIndex"}
      which is currently the same as \code{"numeric"}).}
    \item{\code{rle}:}{Object of class \code{"rle"}, basically a
      \code{\link{list}} with components \code{"lengths"} and
      \code{"values"}, see \code{\link{rle}}.}
  }
}
\section{Methods}{
No methods defined with class "rleDiff" in the signature.
}
\note{
  This is currently an \emph{experimental} auxiliary class
  for the class \code{\linkS4class{abIndex}}, see there.
}
\seealso{
  \code{\link{rle}}, \code{\linkS4class{abIndex}}.
}
\examples{
showClass("rleDiff")
}
\keyword{classes}
