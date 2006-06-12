\name{CHMfactor-class}
\docType{class}
\alias{dCHMsimpl-class}
\alias{lCHMsimpl-class}
\alias{dCHMsuper-class}
\alias{lCHMsuper-class}
\title{CHOLMOD-based Cholesky factorizations}
\description{The hidden class \code{"CHMfactor"} is a class of
  CHOLMOD-based Cholesky factorizations of symmetric, sparse,
  compressed, column-oriented matrices. Such a factorization is
  simplicial (hidden class \code{"CHMsimpl"}) or supernodal (hidden
  class \code{"CHMsuper"}).  Objects that inherit from these classes are either numeric
  factorizations (classes \code{"dCHMsimpl"} and \code{"dCHMsuper"}) or
  symbolic factorizations (classes \code{"dCHMsimpl"} and \code{"dCHMsuper"}).
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("dCHMsuper", ...)}
  but are more commonly created from \code{chol} applied to
  \code{dsCMatrix} or \code{lsCMatrix} objects.
}
\section{Slots}{
  \describe{
    \item{\code{perm}:}{Object of class \code{"integer"} giving the
      permutation of the rows and columns chosen to reduce fill-in and
      for post-ordering.}
    \item{\code{p}:}{Object of class \code{"integer"} of pointers, one
      for each column, to the initial (zero-based) index of elements in
      the column. Only present in classes that contain \code{"CHMsimpl"}.}
    \item{\code{i}:}{Object of class \code{"integer"} of length nnzero
      (number of non-zero elements).  These are the row numbers for
      each non-zero element in the matrix. Only present in classes that
      contain \code{"CHMsimpl"}.}
    \item{\code{x}:}{Object of class \code{"numeric"} - the non-zero
      elements of the matrix.}
  }
}
%\section{Methods}{}
%\references{}
%\author{}
%\note{}
\seealso{
  Class \code{\linkS4class{dgCMatrix}}
}
%\examples{}
\keyword{classes}
\keyword{algebra}
