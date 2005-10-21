\name{symmetricMatrix-class}
\docType{class}
\alias{symmetricMatrix-class}
\title{Virtual Class of Symmetric Matrices in package:Matrix}
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
      inherited from the \code{\link{Matrix-class}}, see there.}
    \item{\code{factors}:}{a list of matrix factorizations, also from the
      \code{Matrix} class.}
  }
}
\section{Extends}{
  Class \code{"Matrix"}, directly.
}
\section{Methods}{
  There's a C function \code{symmetricMatrix_validity()}
  called by the internal validity checking functions.
  %% Currently:
  %% No methods defined with class "symmetricMatrix" in the signature.
}
\seealso{
  \code{\link{triangularMatrix-class}}, and, e.g.,
  \code{\link{dsyMatrix-class}} for numeric \emph{dense} matrices, or
  \code{\link{lsCMatrix-class}} for a logical \emph{sparse} matrix class.
}
\examples{
showClass("symmetricMatrix")

## The names of direct subclasses:
scl <- getClass("symmetricMatrix")@subclasses
directly <- sapply(lapply(scl, slot, "by"), length) == 0
names(scl)[directly]
}
\keyword{classes}
