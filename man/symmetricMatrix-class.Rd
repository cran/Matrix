\name{symmetricMatrix-class}
\docType{class}
\alias{symmetricMatrix-class}
\alias{coerce,matrix,symmetricMatrix-method}
\alias{coerce,denseMatrix,symmetricMatrix-method}
\alias{coerce,CsparseMatrix,symmetricMatrix-method}
\alias{dimnames,symmetricMatrix-method}
%
\title{Virtual Class of Symmetric Matrices in Package Matrix}
\description{
  The virtual class of symmetric matrices, \code{"symmetricMatrix"},
  from the package \pkg{Matrix} contains numeric and logical, dense and
  sparse matrices, e.g., see the examples.

  The main use is in methods (and C functions) that can deal with
  all symmetric matrices.
}
% \section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{uplo}:}{Object of class \code{"character"}.  Must be
      either "U", for upper triangular, and "L", for lower triangular.}
    %% below {Dim, Dimnames} work around Slot parsing buglet (< 2.2.0)
    %% \item{\code{Dim},\code{Dimnames}:}{The dimension (a length-2
    \item{\code{Dim, Dimnames}:}{The dimension (a length-2
      \code{"integer"}) and corresponding names (or \code{NULL}),
      inherited from the \code{\linkS4class{Matrix}}, see there.}
    \item{\code{factors}:}{a list of matrix factorizations, also from the
      \code{Matrix} class.}
  }
}
\section{Extends}{
  Class \code{"Matrix"}, directly.
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "ddiMatrix", to =
	"symmetricMatrix")}: and many other coercion methods, some of
      which are particularly optimized.}
    \item{dimnames}{\code{signature(object = "symmetricMatrix")}:
      returns \emph{symmetric} \code{\link{dimnames}}, even when the
      \code{Dimnames} slot only has row or column names.  This allows to
      save storage for large (typically sparse) symmetric matrices.}
    \item{isSymmetric}{\code{signature(object = "symmetricMatrix")}:
      returns \code{TRUE} trivially.}
  }

  There's a C function \code{symmetricMatrix_validate()}% in ../src/dsyMatrix.c
  called by the internal validity checking functions, and also from
  \code{\link{getValidity}(getClass("symmetricMatrix"))}.
}
\section{Validity and \code{\link{dimnames}}}{
  The validity checks do not require a symmetric \code{Dimnames} slot,
  so it can be \code{list(NULL, <character>)}, e.g., for efficiency.
  However, \code{\link{dimnames}()} and other functions and methods
  should behave as if the dimnames were symmetric, i.e., with both list
  components identical.
}
\seealso{
  \code{\link{isSymmetric}} which has efficient methods
  (\link{isSymmetric-methods}) for the \pkg{Matrix} classes.
  Classes \code{\linkS4class{triangularMatrix}}, and, e.g.,
  \code{\linkS4class{dsyMatrix}} for numeric \emph{dense} matrices, or
  \code{\linkS4class{lsCMatrix}} for a logical \emph{sparse} matrix class.
}
\examples{
showClass("symmetricMatrix")

## The names of direct subclasses:
scl <- getClass("symmetricMatrix")@subclasses
directly <- sapply(lapply(scl, slot, "by"), length) == 0
names(scl)[directly]

## Methods -- applicaple to all subclasses above:
showMethods(classes = "symmetricMatrix")
}
\keyword{classes}
