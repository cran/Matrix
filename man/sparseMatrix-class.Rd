\name{sparseMatrix-class}
\docType{class}
\alias{sparseMatrix-class}
\alias{show,sparseMatrix-method}

\title{Virtual Class "sparseMatrix" --- Mother of Sparse Matrices}
\description{Virtual Mother Class of All Sparse Matrices}
\section{Slots}{
  \code{Dim}, and \code{Dimnames}, see \code{\link{Matrix-class}}.
}
\section{Extends}{
  Class \code{"Matrix"}, directly.
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "sparseMatrix")}: The
      \code{\link{show}} method for sparse matrices prints
      \emph{\dQuote{structural}} zeroes as \code{"."} or blank.}
  }
}
%\author{Martin}
\examples{
showClass("sparseMatrix") ## and look at the help() of its subclasses
}
\keyword{classes}
