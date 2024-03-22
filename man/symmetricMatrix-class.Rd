\name{symmetricMatrix-class}
\title{Virtual Class of Symmetric Matrices in Package Matrix}
%
\docType{class}
\keyword{array}
\keyword{classes}
%
\alias{symmetricMatrix-class}
%
\alias{coerce,matrix,symmetricMatrix-method}
\alias{dimnames,symmetricMatrix-method}
\alias{dimnames<-,symmetricMatrix,NULL-method}
\alias{dimnames<-,symmetricMatrix,list-method}
%
\description{
  The virtual class of symmetric matrices, \code{"symmetricMatrix"},
  from the package \pkg{Matrix} contains numeric and logical, dense and
  sparse matrices, e.g., see the examples with the \dQuote{actual}
  subclasses.

  The main use is in methods (and C functions) that can deal with
  all symmetric matrices, and in \code{as(*, "symmetricMatrix")}.
}
% \section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{Dim, Dimnames}}{inherited from virtual class
      \code{\linkS4class{Matrix}}.  See comments below about
      symmetry of \code{Dimnames}.}
    \item{\code{factors}}{a list of
      \code{\linkS4class{MatrixFactorization}} objects caching
      factorizations of the matrix.  Typically, it is initialized
      as an empty list and updated \dQuote{automagically} whenever
      a factorization is computed.}
    \item{\code{uplo}}{a character string, either \code{"U"} or
      \code{"L"} indicating that only entries in the upper or lower
      triangle are referenced.}
  }
}
\section{Extends}{
  Class \code{"Matrix"}, directly.
}
\section{Methods}{
  \describe{
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
## An example about the symmetric Dimnames:
sy <- sparseMatrix(i= c(2,4,3:5), j= c(4,7:5,5), x = 1:5, dims = c(7,7),
                   symmetric=TRUE, dimnames = list(NULL, letters[1:7]))
sy # shows symmetrical dimnames
sy@Dimnames  # internally only one part is stored
dimnames(sy) # both parts - as sy *is* symmetrical
\dontshow{
local({ nm <- letters[1:7]
  stopifnot(identical(dimnames(sy), list(  nm, nm)),
	    identical(sy@Dimnames , list(NULL, nm)))
})
}%dont
showClass("symmetricMatrix")

## The names of direct subclasses:
scl <- getClass("symmetricMatrix")@subclasses
directly <- sapply(lapply(scl, slot, "by"), length) == 0
names(scl)[directly]

## Methods -- applicaple to all subclasses above:
showMethods(classes = "symmetricMatrix")
}
