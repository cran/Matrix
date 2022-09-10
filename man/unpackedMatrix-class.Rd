\name{unpackedMatrix-class}
\title{Virtual Class \code{"unpackedMatrix"} of Unpacked Dense Matrices}
\Rdversion{1.1}
\docType{class}
%
\alias{unpackedMatrix-class}
%
\alias{coerce,matrix,unpackedMatrix-method}
\alias{coerce,numLike,unpackedMatrix-method}
\alias{diag,unpackedMatrix-method}
\alias{diag<-,unpackedMatrix-method}
\alias{t,unpackedMatrix-method}
%
\description{
  Class \code{"unpackedMatrix"} is the \emph{virtual} class of dense
  matrices in "unpacked" format, storing all \code{m*n} elements of
  an \code{m}-by-\code{n} matrix. It is used to define common methods
  for efficient subsetting, transposing, etc. of its \emph{proper}
  subclasses: currently \code{"[dln]geMatrix"} (unpacked general),
  \code{"[dln]syMatrix"} (unpacked symmetric), \code{"[dln]trMatrix"}
  (unpacked triangular), and subclasses of these, such as
  \code{"\linkS4class{dpoMatrix}"},
  \code{"\linkS4class{Cholesky}"}, and
  \code{"\linkS4class{BunchKaufman}"}.
}
\section{Slots}{
  \describe{
    \item{\code{Dim}, \code{Dimnames}:}{as all \code{\linkS4class{Matrix}} objects.}
  }
}
\section{Extends}{
  Class \code{"\linkS4class{denseMatrix}"}, directly.
  Class \code{"\linkS4class{Matrix}"}, by class \code{"denseMatrix"},
  distance 2.
  Class \code{"mMatrix"}, by class \code{"Matrix"}, distance 3.
  Class \code{"replValueSp"}, by class \code{"Matrix"}, distance 3.
}
\section{Methods}{
  \describe{
    \item{pack}{\code{signature(x = "unpackedMatrix")}: ... }
    \item{unpack}{\code{signature(x = "unpackedMatrix")}: ... }
    \item{isSymmetric}{\code{signature(object = "unpackedMatrix")}: ... }
    \item{isTriangular}{\code{signature(object = "unpackedMatrix")}: ... }
    \item{isDiagonal}{\code{signature(object = "unpackedMatrix")}: ... }
    \item{t}{\code{signature(x = "unpackedMatrix")}: ... }
    \item{diag}{\code{signature(x = "unpackedMatrix")}: ... }
    \item{diag<-}{\code{signature(x = "unpackedMatrix")}: ... }
  }
}
%% \references{
%% }
\author{Mikael Jagan}
%% \note{
%% }
\seealso{
  \code{\link{pack}} and \code{\link{unpack}}; its virtual "complement"
  \code{"\linkS4class{packedMatrix}"}; its proper subclasses
  \code{"\linkS4class{dsyMatrix}"}, \code{"\linkS4class{ltrMatrix}"}, etc.
}
\examples{
showClass("unpackedMatrix")
showMethods(classes = "unpackedMatrix")
}
\keyword{classes}
