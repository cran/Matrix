\name{index-class}
\docType{class}
\alias{index-class}
%% The "[" and "[<-" methods are
%% in ./Xtrct-methods.Rd and ./Subassign-methods.Rd
\title{Virtual Class "index" - Simple Class for Matrix Indices}
\description{The class \code{"index"} is a virtual class used for
  indices (in signatures) for matrix indexing and sub-assignment of
  \pkg{"Matrix"} matrices.

  In fact, it's currently implemented as a simple class union
  (\code{\link{setClassUnion}}) of
  \code{"numeric"},  \code{"logical"} and \code{"character"}.
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
% \section{Methods}{
%   \describe{
%     \item{[<-}{\code{signature(x = "Matrix", i = "missing", j = "missing", value = "index")}: ... }
%     \item{[<-}{\code{signature(x = "denseMatrix", i = "index", j = "missing", value = "numeric")}: ... }
%     \item{[<-}{\code{signature(x = "denseMatrix", i = "index", j = "index", value = "numeric")}: ... }
%     \item{[<-}{\code{signature(x = "denseMatrix", i = "missing", j = "index", value = "numeric")}: ... }
%     \item{[}{\code{signature(x = "Matrix", i = "missing", j = "index", drop = "missing")}: ... }
%     \item{[}{\code{signature(x = "Matrix", i = "index", j = "missing", drop = "missing")}: ... }
%     \item{[}{\code{signature(x = "Matrix", i = "index", j = "index", drop = "missing")}: ... }
%     \item{[}{\code{signature(x = "denseMatrix", i = "index", j = "missing", drop = "logical")}: ... }
%     \item{[}{\code{signature(x = "denseMatrix", i = "index", j = "index", drop = "logical")}: ... }
%     \item{[}{\code{signature(x = "denseMatrix", i = "missing", j = "index", drop = "logical")}: ... }
%     \item{[}{\code{signature(x = "dsparseMatrix", i = "index", j = "missing", drop = "logical")}: ... }
%     \item{[}{\code{signature(x = "dsparseMatrix", i = "index", j = "index", drop = "logical")}: ... }
%     \item{[}{\code{signature(x = "dsparseMatrix", i = "missing", j = "index", drop = "logical")}: ... }
%     \item{[}{\code{signature(x = "lsparseMatrix", i = "index", j = "missing", drop = "logical")}: ... }
%     \item{[}{\code{signature(x = "lsparseMatrix", i = "index", j = "index", drop = "logical")}: ... }
%     \item{[}{\code{signature(x = "lsparseMatrix", i = "missing", j = "index", drop = "logical")}: ... }
%   }
% }
%\author{Martin Maechler}
\seealso{
  \code{\link{[-methods}}, and
  %% FIXME: bug in Rdconv needs '[Matrix]' below:
  \code{\link[Matrix]{Subassign-methods}}, also for examples.
}
\examples{
showClass("index")
}
\keyword{classes}
