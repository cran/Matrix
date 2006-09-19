\name{CHMfactor-class}
\title{CHOLMOD-based Cholesky factorizations}
\docType{class}
\alias{CHMfactor-class}% "mother" virtual
\alias{CHMsimpl-class}% virtual
\alias{CHMsuper-class}% virtual
\alias{dCHMsimpl-class}
\alias{nCHMsimpl-class}
\alias{dCHMsuper-class}
\alias{nCHMsuper-class}
%
\alias{coerce,CHMfactor,sparseMatrix-method}
%
\description{The virtual class \code{"CHMfactor"} is a class of
  CHOLMOD-based Cholesky factorizations of symmetric, sparse,
  compressed, column-oriented matrices.  Such a factorization is
  simplicial (virtual class \code{"CHMsimpl"}) or supernodal (virtual
  class \code{"CHMsuper"}).  Objects that inherit from these classes are
  either numeric factorizations (classes \code{"dCHMsimpl"} and
  \code{"dCHMsuper"}) or symbolic factorizations (classes
  \code{"nCHMsimpl"} and \code{"nCHMsuper"}).
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("dCHMsuper", ...)}
  but are more commonly created via \code{\link{Cholesky}()}.
  applied to  \code{\linkS4class{dsCMatrix}} or
  \code{\linkS4class{lsCMatrix}} objects.
}
\section{Slots}{
  of all \code{"CHMfactor"} and all classes inheriting from it:
  \describe{
    \item{\code{perm}:}{Object of class \code{"integer"} giving the
      permutation of the rows and columns chosen to reduce fill-in and
      for post-ordering.}
    \item{\code{colcount}:}{Object of class \code{"integer"} .... }%% FIXME
    \item{\code{type}:}{Object of class \code{"integer"} .... }
  }

  Slots of the non virtual classes \dQuote{[dl]CHM(super|simpl)}:
  \describe{
    \item{\code{p}:}{Object of class \code{"integer"} of pointers, one
      for each column, to the initial (zero-based) index of elements in
      the column. Only present in classes that contain \code{"CHMsimpl"}.}
    \item{\code{i}:}{Object of class \code{"integer"} of length nnzero
      (number of non-zero elements).  These are the row numbers for
      each non-zero element in the matrix. Only present in classes that
      contain \code{"CHMsimpl"}.}
    \item{\code{x}:}{For the \code{"d*"} classes: \code{"numeric"} - the
      non-zero elements of the matrix.}
  }
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "CHMfactor", to = "sparseMatrix")}}
  }
}
%\references{}
%\author{}
%\note{}
\seealso{
  Class \code{\linkS4class{dgCMatrix}}
}
%\examples{}
\keyword{classes}
\keyword{algebra}
