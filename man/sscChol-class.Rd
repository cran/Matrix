\name{sscChol-class}
\docType{class}
\alias{sscChol-class}
\title{Cholesky decompositions of sscMatrix objects}
\description{The class of Cholesky decompositions of symmetric, sparse
  column matrices from the \code{"sscMatrix"} class.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("sscChol", ...)}
but are more commonly created from \code{chol} applied to an
\code{sscMatrix} object.
}
\section{Slots}{
  \describe{
    \item{\code{perm}:}{Object of class \code{"integer"} giving the
      permutation of the rows and columns chosen to minimize fill-in.
      If pivoting has not been applied this will be an integer vector of
      length 0.}
    \item{\code{diag}:}{A character object indicating if the triangular
      matrix is unit triangular or not.  This slot is inherited from the
      \code{"tscMatrix"} class.  It is always \code{"N"} for this class.}
    \item{\code{uplo}:}{A character object indicating if the upper
      triangle (\code{"U"} or \code{"u"}) or the lower triangle
      (\code{"L"} or \code{"l"}) is stored.  At present only the lower
      triangle form is allowed.}
    \item{\code{p}:}{Object of class \code{"integer"} of pointers, one
      for each column, to the initial (zero-based) index of elements in
      the column.}
    \item{\code{i}:}{Object of class \code{"integer"} of length nnzero
      (number of non-zero elements).  These are the row numbers for
      each non-zero element in the matrix.}
    \item{\code{x}:}{Object of class \code{"numeric"} - the non-zero
      elements of the matrix.}
    \item{\code{factorization}:}{Object of class \code{"list"} - a list
      of factorizations of the matrix. }
    \item{\code{Dim}:}{Object of class \code{"integer"} - the dimensions
     of the matrix - must be an integer vector with exactly two
     non-negative values.}
    \item{\code{D}:}{Object of class \code{"numeric"} - the diagonal
      elements of the matrix D in the LDL' decomposition.}
    \item{\code{Parent}:}{Object of class \code{"integer"} - the
      elimination tree of the (possibly permuted) symmetric matrix.}
  }
}
\section{Extends}{
Class \code{"tscMatrix"}, directly.
Class \code{"cscMatrix"}, by class \code{"tscMatrix"}.
}
\section{Methods}{
No methods defined with class "sscChol" in the signature.
}
%\references{}
%\author{}
%\note{}
\seealso{
  \code{\link{cscMatrix-class}}
}
%\examples{}
\examples{
data(mm)
xpx = crossprod(mm)
str(xpx)
str(chol(xpx))
}
\keyword{classes}
\keyword{algebra}
