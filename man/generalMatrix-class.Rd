\name{generalMatrix-class}
\title{Class "generalMatrix" of General Matrices}
%
\docType{class}
\keyword{array}
\keyword{classes}
%
\alias{generalMatrix-class}
%
\alias{coerce,generalMatrix,packedMatrix-method}
\alias{coerce,matrix,generalMatrix-method}
\alias{coerce,vector,generalMatrix-method}
\alias{dimnames<-,generalMatrix,NULL-method}
\alias{dimnames<-,generalMatrix,list-method}
%
\description{
  Virtual class of \dQuote{general} matrices; i.e., matrices
  that do not have a known property such as symmetric, triangular, or
  diagonal.
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{Dim, Dimnames}}{inherited from virtual class
      \code{\linkS4class{Matrix}}.}
    \item{\code{factors}}{a list of
      \code{\linkS4class{MatrixFactorization}} objects caching
      factorizations of the matrix.  Typically, it is initialized
      as an empty list and updated \dQuote{automagically} whenever
      a factorization is computed.}
  }
}
\section{Extends}{
  Class \code{"Matrix"}, directly.
}
% \section{Methods}{
% No methods defined with class "generalMatrix" in the signature.
% }
\seealso{
  Virtual classes
  \code{\linkS4class{symmetricMatrix}},
  \code{\linkS4class{triangularMatrix}}, and
  \code{\linkS4class{diagonalMatrix}}.
}
% \examples{
% ##---- Should be DIRECTLY executable !! ----
% }
