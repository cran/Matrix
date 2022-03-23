\name{packedMatrix-class}
\title{Class \code{"packedMatrix"} of packed dense symmetric/triangular Matrices}
\Rdversion{1.1}
\docType{class}
\alias{packedMatrix-class}
%
\alias{[,packedMatrix,missing,missing,missing-method}
\alias{[,packedMatrix,NULL,missing,missing-method}
\alias{[,packedMatrix,index,missing,missing-method}
\alias{[,packedMatrix,matrix,missing,missing-method}
\alias{[,packedMatrix,lMatrix,missing,missing-method}
\alias{[,packedMatrix,missing,NULL,missing-method}
\alias{[,packedMatrix,NULL,NULL,missing-method}
\alias{[,packedMatrix,index,NULL,missing-method}
\alias{[,packedMatrix,matrix,NULL,missing-method}
\alias{[,packedMatrix,lMatrix,NULL,missing-method}
\alias{[,packedMatrix,missing,index,missing-method}
\alias{[,packedMatrix,NULL,index,missing-method}
\alias{[,packedMatrix,index,index,missing-method}
\alias{[,packedMatrix,matrix,index,missing-method}
\alias{[,packedMatrix,lMatrix,index,missing-method}
\alias{[,packedMatrix,missing,missing,logical-method}
\alias{[,packedMatrix,NULL,missing,logical-method}
\alias{[,packedMatrix,index,missing,logical-method}
\alias{[,packedMatrix,matrix,missing,logical-method}
\alias{[,packedMatrix,lMatrix,missing,logical-method}
\alias{[,packedMatrix,missing,NULL,logical-method}
\alias{[,packedMatrix,NULL,NULL,logical-method}
\alias{[,packedMatrix,index,NULL,logical-method}
\alias{[,packedMatrix,matrix,NULL,logical-method}
\alias{[,packedMatrix,lMatrix,NULL,logical-method}
\alias{[,packedMatrix,missing,index,logical-method}
\alias{[,packedMatrix,NULL,index,logical-method}
\alias{[,packedMatrix,index,index,logical-method}
\alias{[,packedMatrix,matrix,index,logical-method}
\alias{[,packedMatrix,lMatrix,index,logical-method}
%
\alias{diag,packedMatrix-method}
\alias{diag<-,packedMatrix-method}
%
\alias{t,packedMatrix-method}
\description{
  Class \code{"packedMatrix"} is a \emph{virtual} class of packed dense
  symmetric or triangular Matrices.

  It is used to define common methods for efficient subsetting, transposing
  etc for the \emph{actual} classes it contains, i.e.,  the currently 2 x 3
  classes \code{"[dln]spMatrix"} (packed symmetric) and
  \code{"[dln]tpMatrix"} (packed triangular), and subclasses of these, such as
  \code{"\linkS4class{dppMatrix}"},
  \code{"\linkS4class{pCholesky}"}, and \code{"\linkS4class{pBunchKaufman}"}.
}
\section{Slots}{
  \describe{
    \item{\code{uplo}:}{\code{"character"}; either "U", for upper triangular, and "L", for lower.}
    \item{\code{Dim}, \code{Dimnames}:}{as all \code{\linkS4class{Matrix}} objects.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{Matrix}"}, directly.
Class \code{"mMatrix"}, by class "Matrix", distance 2.
Class \code{"\linkS4class{replValue}"}, by class "Matrix", distance 2.
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "packedMatrix", i = "missing", j = "missing", drop = "missing")}: ... }
    \item{[}{\code{signature(x = "packedMatrix", i = "NULL", j = "missing", drop = "missing")}: ... }
    \item{[}{\code{signature(x = "packedMatrix", i = "index", j = "missing", drop = "missing")}: ... }
    \item{[}{\code{signature(x = "packedMatrix", i = "matrix", j = "missing", drop = "missing")}: ... }
    \item{[}{\code{signature(x = "packedMatrix", i = "lMatrix", j = "missing", drop = "missing")}: ... }
    \item{[}{\code{signature(x = "packedMatrix", i = "missing", j = "NULL", drop = "missing")}: ... }
    \item{[}{\code{signature(x = "packedMatrix", i = "NULL", j = "NULL", drop = "missing")}: ... }
    \item{[}{\code{signature(x = "packedMatrix", i = "index", j = "NULL", drop = "missing")}: ... }
    \item{[}{\code{signature(x = "packedMatrix", i = "matrix", j = "NULL", drop = "missing")}: ... }
    \item{[}{\code{signature(x = "packedMatrix", i = "lMatrix", j = "NULL", drop = "missing")}: ... }
    \item{[}{\code{signature(x = "packedMatrix", i = "missing", j = "index", drop = "missing")}: ... }
    \item{[}{\code{signature(x = "packedMatrix", i = "NULL", j = "index", drop = "missing")}: ... }
    \item{[}{\code{signature(x = "packedMatrix", i = "index", j = "index", drop = "missing")}: ... }
    \item{[}{\code{signature(x = "packedMatrix", i = "matrix", j = "index", drop = "missing")}: ... }
    \item{[}{\code{signature(x = "packedMatrix", i = "lMatrix", j = "index", drop = "missing")}: ... }
    \item{[}{\code{signature(x = "packedMatrix", i = "missing", j = "missing", drop = "logical")}: ... }
    \item{[}{\code{signature(x = "packedMatrix", i = "NULL", j = "missing", drop = "logical")}: ... }
    \item{[}{\code{signature(x = "packedMatrix", i = "index", j = "missing", drop = "logical")}: ... }
    \item{[}{\code{signature(x = "packedMatrix", i = "matrix", j = "missing", drop = "logical")}: ... }
    \item{[}{\code{signature(x = "packedMatrix", i = "lMatrix", j = "missing", drop = "logical")}: ... }
    \item{[}{\code{signature(x = "packedMatrix", i = "missing", j = "NULL", drop = "logical")}: ... }
    \item{[}{\code{signature(x = "packedMatrix", i = "NULL", j = "NULL", drop = "logical")}: ... }
    \item{[}{\code{signature(x = "packedMatrix", i = "index", j = "NULL", drop = "logical")}: ... }
    \item{[}{\code{signature(x = "packedMatrix", i = "matrix", j = "NULL", drop = "logical")}: ... }
    \item{[}{\code{signature(x = "packedMatrix", i = "lMatrix", j = "NULL", drop = "logical")}: ... }
    \item{[}{\code{signature(x = "packedMatrix", i = "missing", j = "index", drop = "logical")}: ... }
    \item{[}{\code{signature(x = "packedMatrix", i = "NULL", j = "index", drop = "logical")}: ... }
    \item{[}{\code{signature(x = "packedMatrix", i = "index", j = "index", drop = "logical")}: ... }
    \item{[}{\code{signature(x = "packedMatrix", i = "matrix", j = "index", drop = "logical")}: ... }
    \item{[}{\code{signature(x = "packedMatrix", i = "lMatrix", j = "index", drop = "logical")}: ... }

    \item{diag}{\code{signature(x = "packedMatrix")}: ... }
    \item{diag<-}{\code{signature(x = "packedMatrix")}: ... }

    \item{t}{\code{signature(x = "packedMatrix")}: ... }
  }
}
%% \references{
%% }
\author{Mikael Jagan}
%% \note{
%% }
\seealso{
  \code{\link{pack}} and \code{\link{unpack}};
  its subclasses, \code{"\linkS4class{dtpMatrix}"},
  \code{"\linkS4class{lspMatrix}"}, etc.
}
\examples{
showClass("packedMatrix")
showMethods(class = "packedMatrix")
}
\keyword{classes}
