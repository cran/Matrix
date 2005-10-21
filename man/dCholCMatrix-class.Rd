\name{dCholCMatrix-class}
\docType{class}
\alias{dCholCMatrix-class}
\alias{lCholCMatrix-class}
\alias{solve,lCholCMatrix,missing-method}
\alias{solve,lCholCMatrix,lgCMatrix-method}
\alias{t,lCholCMatrix-method}
\title{Cholesky Decompositions of dsCMatrix Objects}
\description{\code{"dCholCMatrix"} is class of Cholesky decompositions
  of symmetric, sparse, compressed, column-oriented matrices (the
  \code{"dsCMatrix"} class).  \code{"lCholCMatrix"} is class showing the
  logical structure (positions of the non-zeros) of these Cholesky decompositions.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("dCholCMatrix", ...)}
  but are more commonly created from \code{chol} applied to
  \code{dsCMatrix} or \code{lsCMatrix} objects.
}
\section{Slots}{
  \describe{
    \item{\code{perm}:}{Object of class \code{"integer"} giving the
      permutation of the rows and columns chosen to minimize fill-in.
      If pivoting has not been applied this will be an integer vector of
      length 0.}
    \item{\code{diag}:}{A character object indicating if the triangular
      matrix is unit triangular or not.  This slot is inherited from the
      \code{"dtCMatrix"} class.  It is always \code{"N"} for this class.}
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
      elements of the matrix.  Present only in the \code{"dCholCMatrix"}
      class.}
    \item{\code{factors}:}{Object of class \code{"list"} - a list
      of factorizations of the matrix. }
    \item{\code{Dim}:}{Object of class \code{"integer"} - the dimensions
     of the matrix - must be an integer vector with exactly two
     non-negative values.}
    \item{\code{Dimnames}:}{list of length two; inherited from the
      \code{"Matrix"} class, see \code{\link{Matrix-class}}.}
    \item{\code{D}:}{Object of class \code{"numeric"} - the diagonal
      elements of the matrix D in the LDL' decomposition. Present only
      in the \code{"dCholCMatrix"} class.}
    \item{\code{Parent}:}{Object of class \code{"integer"} - the
      elimination tree of the (possibly permuted) symmetric matrix.}
  }
}
\section{Extends}{
  Class \code{"dtCMatrix"}, directly.
  Class \code{"dgCMatrix"}, by class \code{"dtCMatrix"}.
}
\section{Methods}{
  No methods defined with class "dCholCMatrix" in the signature.
}
%\references{}
%\author{}
%\note{}
\seealso{
  \code{\link{dgCMatrix-class}}
}
\examples{
data(mm)
xpx <- crossprod(mm)
str(xpx)
str(ch <- chol(xpx))

## Visualize the sparseness:
dq <- function(ch) paste('"',ch,'"', sep="") ## dQuote(<UTF-8>) gives bad plots
image(xpx, main=paste("crossprod(mm) : Sparse", dq(class(xpx))))
image(ch, main= paste("chol (crossprod(mm)) : Sparse", dq(class(ch))))
}
\keyword{classes}
\keyword{algebra}
