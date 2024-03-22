\name{packedMatrix-class}
\title{Virtual Class \code{"packedMatrix"} of Packed Dense Matrices}
%
\docType{class}
\keyword{array}
\keyword{classes}
%
\alias{packedMatrix-class}
%
\alias{coerce,matrix,packedMatrix-method}
\alias{cov2cor,packedMatrix-method}
%
\description{
  Class \code{"packedMatrix"} is the \emph{virtual} class of dense
  symmetric or triangular matrices in "packed" format, storing only
  the \code{choose(n+1,2) == n*(n+1)/2} elements of the upper or
  lower triangle of an \code{n}-by-\code{n} matrix. It is used to
  define common methods for efficient subsetting, transposing, etc.
  of its \emph{proper} subclasses: currently \code{"[dln]spMatrix"}
  (packed symmetric), \code{"[dln]tpMatrix"} (packed triangular),
  and subclasses of these, such as
  \code{"\linkS4class{dppMatrix}"}.
}
\section{Slots}{
  \describe{
    \item{\code{uplo}:}{\code{"character"}; either "U", for upper triangular, and "L", for lower.}
    \item{\code{Dim}, \code{Dimnames}:}{as all \code{\linkS4class{Matrix}} objects.}
  }
}
\section{Extends}{
  Class \code{"\linkS4class{denseMatrix}"}, directly.
  Class \code{"\linkS4class{Matrix}"}, by class \code{"denseMatrix"},
  distance 2.
}
\section{Methods}{
  \describe{
    \item{pack}{\code{signature(x = "packedMatrix")}: ... }
    \item{unpack}{\code{signature(x = "packedMatrix")}: ... }
    \item{isSymmetric}{\code{signature(object = "packedMatrix")}: ... }
    \item{isTriangular}{\code{signature(object = "packedMatrix")}: ... }
    \item{isDiagonal}{\code{signature(object = "packedMatrix")}: ... }
    \item{t}{\code{signature(x = "packedMatrix")}: ... }
    \item{diag}{\code{signature(x = "packedMatrix")}: ... }
    \item{diag<-}{\code{signature(x = "packedMatrix")}: ... }
  }
}
%% \references{
%% }
\author{Mikael Jagan}
%% \note{
%% }
\seealso{
  \code{\link{pack}} and \code{\link{unpack}}; its virtual "complement"
  \code{"\linkS4class{unpackedMatrix}"}; its proper subclasses
  \code{"\linkS4class{dspMatrix}"}, \code{"\linkS4class{ltpMatrix}"}, etc.
}
\examples{
showClass("packedMatrix")
showMethods(classes = "packedMatrix")
}
