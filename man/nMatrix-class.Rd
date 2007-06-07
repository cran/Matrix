\name{nMatrix-class}
\docType{class}
\alias{nMatrix-class}
\alias{show,nMatrix-method}
\alias{coerce,nMatrix,matrix-method}
\alias{coerce,nMatrix,dMatrix-method}
\alias{coerce,dMatrix,nMatrix-method}
\alias{coerce,nMatrix,lMatrix-method}
\alias{coerce,lMatrix,nMatrix-method}
% %\alias{coerce,dMatrix,dgeMatrix-method}
% \alias{[,dMatrix,lMatrix,missing,ANY-method}
% \alias{[,dMatrix,logical,missing,ANY-method}
% % Group methods
% \alias{Arith,dMatrix,dMatrix-method}
% \alias{Math2,dMatrix,numeric-method}% Math2 = round + signif, but
% \alias{Math2,dMatrix,missing-method}
% \alias{Summary,dMatrix-method}
% \alias{Compare,dMatrix,dMatrix-method}
% \alias{Compare,dMatrix,numeric-method}
% \alias{Compare,numeric,dMatrix-method}
% % for silly reasons, need these 2+3 as well:
% \alias{round,dMatrix,numeric-method}
% \alias{signif,dMatrix,numeric-method}
% \alias{log,dMatrix-method}
% \alias{gamma,dMatrix-method}
% \alias{lgamma,dMatrix-method}
% %
% \alias{zapsmall,dMatrix-method}
%
\title{Class "nMatrix" of Non-zero Pattern Matrices}
\description{
  The \code{lMatrix} class is the virtual \dQuote{mother} class of all
  \emph{\bold{n}on-zero pattern} matrices in the \pkg{Matrix} package.
}
%\section{Objects from the Class}{A virtual Class: No objects may be
%  created from it.
%}
\section{Slots}{
  Common to \emph{all} matrix object in the package:
  \describe{
    \item{\code{Dim}:}{Object of class \code{"integer"} - the dimensions
     of the matrix - must be an integer vector with exactly two
     non-negative values.}
    \item{\code{Dimnames}:}{list of length two; each component
      containing NULL or a \code{\link{character}} vector length
      equal the corresponding \code{Dim} element.}
  }
}
\section{Methods}{
%   There are (relatively simple) group methods (see, e.g., \code{\link{Arith}})
%   \describe{
%     \item{Arith}{\code{signature(e1 = "dMatrix", e2 = "dMatrix")}: ... }
%     \item{Arith}{\code{signature(e1 = "dMatrix", e2 = "numeric")}: ... }
%     \item{Arith}{\code{signature(e1 = "numeric", e2 = "dMatrix")}: ... }
%     \item{Math}{\code{signature(x = "dMatrix")}: ... }
%     \item{Math2}{\code{signature(x = "dMatrix", digits = "numeric")}:
%       this group contains \code{\link{round}()} and \code{\link{signif}()}.}
%     \item{Compare}{\code{signature(e1 = "numeric", e2 = "dMatrix")}: ... }
%     \item{Compare}{\code{signature(e1 = "dMatrix", e2 = "numeric")}: ... }
%     \item{Compare}{\code{signature(e1 = "dMatrix", e2 = "dMatrix")}: ... }
%     \item{Summary}{\code{signature(x = "dMatrix")}: The \code{"Summary"}
%       group contains the seven functions
%       \code{\link{max}()}, \code{\link{min}()}, \code{\link{range}()},
%       \code{\link{prod}()}, \code{\link{sum}()},
%       \code{\link{any}()}, and \code{\link{all}()}.}
%   }
% The following methods are defined for all n* matrices:
  \describe{
    \item{coerce}{\code{signature(from = "nMatrix", to = "matrix")}: ... }
  }
}
%\references{}
%\author{Martin + Doug}
\seealso{
%   The classes \code{\linkS4class{dgeMatrix}},
%   \code{\linkS4class{dgCMatrix}}, and
%% FIXME
\code{\linkS4class{Matrix}}.
}
\examples{
 showClass("nMatrix")

}
\keyword{classes}
\keyword{algebra}
