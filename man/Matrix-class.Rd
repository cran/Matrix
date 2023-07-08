\name{Matrix-class}
\title{Virtual Class "Matrix" of Matrices}
%
\docType{class}
\keyword{array}
\keyword{classes}
%
\alias{Matrix-class}
%
\alias{!,Matrix-method}
\alias{&,Matrix,ddiMatrix-method}
\alias{&,Matrix,ldiMatrix-method}
\alias{*,Matrix,ddiMatrix-method}
\alias{*,Matrix,ldiMatrix-method}
\alias{+,Matrix,missing-method}
\alias{-,Matrix,missing-method}
\alias{Arith,Matrix,Matrix-method}
\alias{Arith,Matrix,lsparseMatrix-method}
\alias{Arith,Matrix,nsparseMatrix-method}
\alias{Logic,ANY,Matrix-method}
\alias{Logic,Matrix,ANY-method}
\alias{Logic,Matrix,nMatrix-method}
\alias{Math2,Matrix-method}
\alias{Ops,ANY,Matrix-method}
\alias{Ops,Matrix,ANY-method}
\alias{Ops,Matrix,NULL-method}
\alias{Ops,Matrix,ddiMatrix-method}
\alias{Ops,Matrix,ldiMatrix-method}
\alias{Ops,Matrix,matrix-method}
\alias{Ops,Matrix,sparseVector-method}
\alias{Ops,NULL,Matrix-method}
\alias{Ops,matrix,Matrix-method}
\alias{Summary,Matrix-method}
\alias{^,Matrix,ddiMatrix-method}
\alias{^,Matrix,ldiMatrix-method}
\alias{as.array,Matrix-method}
\alias{as.logical,Matrix-method}
\alias{as.matrix,Matrix-method}
\alias{as.numeric,Matrix-method}
\alias{as.vector,Matrix-method}
\alias{coerce,ANY,Matrix-method}
\alias{coerce,Matrix,complex-method}
\alias{coerce,Matrix,corMatrix-method}
\alias{coerce,Matrix,diagonalMatrix-method}
\alias{coerce,Matrix,dpoMatrix-method}
\alias{coerce,Matrix,dppMatrix-method}
\alias{coerce,Matrix,indMatrix-method}
\alias{coerce,Matrix,integer-method}
\alias{coerce,Matrix,logical-method}
\alias{coerce,Matrix,numeric-method}
\alias{coerce,Matrix,pMatrix-method}
\alias{coerce,Matrix,pcorMatrix-method}
\alias{coerce,Matrix,symmetricMatrix-method}
\alias{coerce,Matrix,triangularMatrix-method}
\alias{coerce,Matrix,vector-method}
\alias{coerce,matrix,Matrix-method}
\alias{determinant,Matrix,missing-method}
\alias{determinant,Matrix,logical-method}
\alias{diff,Matrix-method}
\alias{dim,Matrix-method}
\alias{dimnames,Matrix-method}
\alias{dimnames<-,Matrix,NULL-method}
\alias{dimnames<-,Matrix,list-method}
\alias{drop,Matrix-method}
\alias{head,Matrix-method}
\alias{initialize,Matrix-method}
\alias{length,Matrix-method}
\alias{tail,Matrix-method}
\alias{unname,Matrix-method}
%
\alias{det} % for our copy of base::det
%
\description{
  The \code{Matrix} class is a class contained by all actual
  classes in the \pkg{Matrix} package.  It is a \dQuote{virtual} class.
}
\section{Slots}{
  \describe{
    \item{\code{Dim}}{an integer vector of length 2 giving the
      dimensions of the matrix.}
    \item{\code{Dimnames}}{a list of length 2.  Each element must
      be \code{NULL} or a character vector of length equal to the
      corresponding element of \code{Dim}.}
  }
}
\section{Methods}{
  \describe{
    \item{determinant}{\code{signature(x = "Matrix", logarithm = "missing")}: and}
    \item{determinant}{\code{signature(x = "Matrix", logarithm = "logical")}:
      compute the (\eqn{\log}) determinant of \code{x}.  The method
      chosen depends on the actual Matrix class of \code{x}.  Note that
      \code{\link[base]{det}} also works for all our matrices, calling the
      appropriate \code{determinant()} method.  The \code{Matrix::det}
      is an exact copy of \code{base::det}, but in the correct
      namespace, and hence calling the S4-aware version of
      \code{determinant()}.).}
    \item{diff}{\code{signature(x = "Matrix")}: As \code{\link{diff}()}
      for traditional matrices, i.e., applying \code{diff()} to each
      column.}
    \item{dim}{\code{signature(x = "Matrix")}: extract matrix dimensions
      \code{\link{dim}}.}
    \item{dim<-}{\code{signature(x = "Matrix", value = "ANY")}: where
      \code{value} is integer of length 2.  Allows to \emph{reshape}
      Matrix objects, but only when \code{prod(value) == prod(dim(x))}.}
    \item{dimnames}{\code{signature(x = "Matrix")}: extract
      \code{\link{dimnames}}.}
    \item{dimnames<-}{\code{signature(x = "Matrix", value = "list")}: set
      the \code{dimnames} to a \code{\link{list}} of length 2, see
      \code{\link{dimnames<-}}.}
    \item{length}{\code{signature(x = "Matrix")}: simply defined as
      \code{prod(dim(x))} (and hence of mode \code{"double"}).}
    \item{show}{\code{signature(object = "Matrix")}: \code{\link{show}}
      method for \code{\link{print}}ing.  For printing \emph{sparse}
      matrices, see \code{\link{printSpMatrix}}.}
    \item{image}{\code{signature(object = "Matrix")}: draws an
      \code{\link{image}} of the matrix entries, using
      \code{\link[lattice]{levelplot}()} from package \pkg{lattice}.}
    \item{head}{\code{signature(object = "Matrix")}: return only the
      \emph{\dQuote{head}}, i.e., the first few rows.}
    \item{tail}{\code{signature(object = "Matrix")}: return only the
      \emph{\dQuote{tail}}, i.e., the last few rows of the respective matrix.}
    \cr %------------------------------------
    \item{as.matrix, as.array}{\code{signature(x = "Matrix")}: the same as
      \code{as(x, "matrix")}; see also the note below.}
    \item{as.vector}{\code{signature(x = "Matrix", mode = "missing")}:
      \code{as.vector(m)} should be identical to \code{as.vector(as(m,
	"matrix"))}, implemented more efficiently for some subclasses.}
    \item{as(x, "vector"), as(x, "numeric")}{etc, similarly.}
    \item{coerce}{\code{signature(from = "ANY", to = "Matrix")}: This
      relies on a correct \code{\link{as.matrix}()} method for \code{from}.}
  }
  There are many more methods that (conceptually should) work for all
  \code{"Matrix"} objects, e.g., \code{\link{colSums}},
  \code{\link{rowMeans}}.  Even \pkg{base} functions may work
  automagically (if they first call \code{\link{as.matrix}()} on their
  principal argument), e.g., \code{\link{apply}}, \code{\link{eigen}},
  \code{\link{svd}} or \code{\link{kappa}} all do work via coercion to a
  \dQuote{traditional} (dense) \code{\link{matrix}}.
  %% --> ../tests/base-matrix-fun.R
}
\note{
  Loading the \code{Matrix} namespace \dQuote{overloads}
  \code{\link{as.matrix}} and \code{\link{as.array}} in the \pkg{base}
  namespace by the equivalent of \code{function(x) as(x, "matrix")}.
  Consequently, \code{as.matrix(m)} or \code{as.array(m)} will properly
  work when \code{m} inherits from the \code{"Matrix"} class ---
  \emph{also} for functions in package \pkg{base} and other packages.
  E.g., \code{\link{apply}} or \code{\link{outer}} can therefore be applied
  to \code{"Matrix"} matrices.
}
\author{Douglas Bates \email{bates@stat.wisc.edu} and Martin Maechler}
\seealso{
  the classes \code{\linkS4class{dgeMatrix}},
  \code{\linkS4class{dgCMatrix}}, and function
  \code{\link{Matrix}} for construction (and examples).

  Methods, e.g., for \code{\link[=kronecker-methods]{kronecker}}.
}
\examples{
slotNames("Matrix")

cl <- getClass("Matrix")
names(cl@subclasses) # more than 40 ..

showClass("Matrix")#> output with slots and all subclasses

(M <- Matrix(c(0,1,0,0), 6, 4))
dim(M)
diag(M)
cm <- M[1:4,] + 10*Diagonal(4)
diff(M)
## can reshape it even :
dim(M) <- c(2, 12)
M
stopifnot(identical(M, Matrix(c(0,1,0,0), 2,12)),
          all.equal(det(cm),
                    determinant(as(cm,"matrix"), log=FALSE)$modulus,
                    check.attributes=FALSE))
}
