\name{triangularMatrix-class}
\docType{class}
\alias{triangularMatrix-class}
\title{Virtual Class of Triangular Matrices in package:Matrix}
\description{
  The virtual class of triangular matrices,\code{"triangularMatrix"},
  from the package \pkg{Matrix} contains numeric and logical, dense and
  sparse matrices, e.g., see the examples.

  The main use will be in methods (and C functions) that can deal with
  all triangular matrices.
}
% \section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{uplo}:}{Object of class \code{"character"}.  Must be
      either "U", for upper triangular, and "L", for lower triangular.}
    \item{\code{diag}:}{Object of class \code{"character"}.  Must be
      either \code{"U"}, for unit triangular (diagonal is all ones), or
      \code{"N"} for non-unit.  The diagonal elements are not
      accessed internally when \code{diag} is \code{"U"} (but need to be
      allocated, i.e., \code{length(x)} does not depend on \code{diag}).}
    %% below {Dim, Dimnames} work around Slot parsing buglet (< 2.2.0)
    \item{\code{Dim, Dimnames}:}{The dimension (a length-2
      \code{"integer"}) and corresponding names (or \code{NULL}),
      inherited from the \code{\link{Matrix-class}}, see there.}
    \item{\code{factors}:}{a list of matrix factorizations, also from the
      \code{Matrix} class.}
  }
}
\section{Extends}{
  Class \code{"Matrix"}, directly.
}
\section{Methods}{
  There's a C function \code{triangularMatrix_validity()}
  called by the internal validity checking functions.
  %% Currently:
  %% No methods defined with class "triangularMatrix" in the signature.
}
\seealso{
  \code{\link{symmetricMatrix-class}}, and, e.g.,
  \code{\link{dtrMatrix-class}} for numeric \emph{dense} matrices, or
  \code{\link{ltCMatrix-class}} for a logical \emph{sparse} matrix
  subclass of \code{"triangularMatrix"}.
}
\examples{
showClass("triangularMatrix")

## The names of direct subclasses:
scl <- getClass("triangularMatrix")@subclasses
directly <- sapply(lapply(scl, slot, "by"), length) == 0
names(scl)[directly]
}
\keyword{classes}
