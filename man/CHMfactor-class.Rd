\name{CHMfactor-class}
\title{CHOLMOD-based Cholesky Factorizations}
\docType{class}
\alias{CHMfactor-class}% "mother" virtual
\alias{CHMsimpl-class}% virtual
\alias{CHMsuper-class}% virtual
\alias{dCHMsimpl-class}
\alias{nCHMsimpl-class}
\alias{dCHMsuper-class}
\alias{nCHMsuper-class}
%
\alias{coerce,CHMfactor,Matrix-method}
\alias{coerce,CHMfactor,sparseMatrix-method}
\alias{solve,CHMfactor,ddenseMatrix-method}
\alias{solve,CHMfactor,dsparseMatrix-method}
\alias{solve,CHMfactor,numeric-method}
\alias{solve,CHMfactor,matrix-method}
\alias{solve,CHMfactor,ANY-method}
\alias{determinant,CHMfactor,missing-method}
\alias{determinant,CHMfactor,logical-method}
\alias{update,CHMfactor-method}
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
\usage{% want usage for the update method which has "surprising arguments"
\S4method{update}{CHMfactor}(object, parent, mult = 0, \dots)
}
\arguments{
  \item{object}{a \code{"CHMfactor"} object.}
  \item{parent}{a \code{"\linkS4class{dsCMatrix}"} matrix or
    \code{"\linkS4class{dgCMatrix}"} object with the same nonzero
    pattern as the matrix that generated \code{object}.  If
    \code{parent} is a symmetric matrix (class
    \code{"\linkS4class{dsCMatrix}"}) then \code{object} should be a
    decomposition of a matrix with the same nonzero pattern as
    \code{parent}.  If \code{parent} is not symmetric then \code{object}
    should be the decomposition of a matrix with the same nonzero
    pattern as \code{tcrossprod(parent)}.}
  \item{mult}{a numeric scalar (default 0).  \code{mult} times the
    identity matrix is (implicitly) added to \code{parent} or
    \code{tcrossprod(parent)} before updating the decomposition
    \code{object}.}
  \item{\dots}{Optional additional arguments to the methods.}
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("dCHMsuper", ...)}
  but are more commonly created via \code{\link{Cholesky}()},
  applied to  \code{\linkS4class{dsCMatrix}} or
  \code{\linkS4class{lsCMatrix}} objects.
}
\section{Slots}{
  of \code{"CHMfactor"} and all classes inheriting from it:
  \describe{
    \item{\code{perm}:}{An integer vector giving the 0-based
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
      each non-zero element in the matrix.  Only present in classes that
      contain \code{"CHMsimpl"}.}
    \item{\code{x}:}{For the \code{"d*"} classes: \code{"numeric"} - the
      non-zero elements of the matrix.}
  }
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "CHMfactor", to = "sparseMatrix")}}
    \item{image}{\code{signature(x = "CHMfactor")}}

    \item{solve}{\code{signature(a = "CHMfactor", b = "ddenseMatrix")}
      The \code{solve} methods for a \code{"CHMfactor"} object take an
      optional third argument \code{system} whose value can be one of the
      character strings \code{"A"}, \code{"LDLt"}, \code{"LD"},
      \code{"DLt"}, \code{"L"}, \code{"Lt"}, \code{"D"}, \code{"P"} or
      \code{"Pt"}.  This argument describes the system to be solved.  The
      default is to solve \code{Ax = b} for \code{x} where \code{A} is
      sparse, positive-definite matrix that was factored to produce
      \code{"a"}.}
    \item{solve}{\code{signature(a = "CHMfactor", b = "matrix")}, and}
    \item{solve}{\code{signature(a = "CHMfactor", b = "numeric")} call
      the above \code{signature(*, b = "ddenseMatrix")} method.}

    \item{solve}{\code{signature(a = "CHMfactor", b = "dsparseMatrix")}
      takes the same third argument \code{system} and calls the
      corresponding sparse CHOLMOD algorithm.}

    \item{determinant}{\code{signature(x = "CHMfactor", logarithm =
      "logical")} returns the determinant (or the logarithm of the
      determinant, if \code{logarithm = TRUE}, the default) of \emph{the
      matrix whose decomposition is represented by} \code{x}.  The
      determinant of the decomposition is not well defined because the
      decomposition may be \dQuote{LDL'} or \dQuote{LL'}.  The
      determinant of the original matrix that was decomposed is well
      defined.}

    \item{update}{\code{signature(object = "CHMfactor")}.  The
      \code{\link{update}} method requires an additional argument
      \code{parent}, which is a \code{"\linkS4class{dsCMatrix}"} object with
      the same structure of nonzeros as the matrix that was decomposed
      to produce \code{object}, and provides an optional argument
      \code{mult}, a numeric scalar.  This method updates the numeric
      values in \code{object} to the decomposition of \eqn{A+mI} where
      \eqn{A} is the matrix represented by \code{parent} and \eqn{m} is
      the scalar \code{mult}. Because only the numeric values are
      updated this method should be faster than creating and decomposing
      \eqn{A+mI}.  It is not uncommon to want, say, the determinant of
      \eqn{A+mI} for many different values of \eqn{m}.  This method
      would be the preferred approach in such cases.}

  }
}
%\references{}
%\author{}
%\note{}
\seealso{
  \code{\link{Cholesky}}, also for examples;
  class \code{\linkS4class{dgCMatrix}}.
}
%\examples{} are in ./Cholesky.Rd
\keyword{classes}
\keyword{algebra}
