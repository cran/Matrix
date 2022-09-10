\name{lsparseMatrix-classes}
\title{Sparse logical matrices}
\docType{class}
%
\alias{lsparseMatrix-class}
\alias{lgCMatrix-class}
\alias{lgRMatrix-class}
\alias{lgTMatrix-class}
\alias{ltCMatrix-class}
\alias{ltRMatrix-class}
\alias{ltTMatrix-class}
\alias{lsCMatrix-class}
\alias{lsRMatrix-class}
\alias{lsTMatrix-class}
% lsparse
\alias{!,lsparseMatrix-method}
\alias{Arith,lsparseMatrix,Matrix-method}
\alias{Logic,lsparseMatrix,ldenseMatrix-method}
\alias{Logic,lsparseMatrix,lsparseMatrix-method}
\alias{Ops,lsparseMatrix,lsparseMatrix-method}
\alias{Ops,lsparseMatrix,nsparseMatrix-method}
\alias{as.logical,lsparseMatrix-method}
\alias{as.numeric,lsparseMatrix-method}
\alias{coerce,lsparseMatrix,dMatrix-method}
\alias{coerce,lsparseMatrix,dsparseMatrix-method}
\alias{coerce,lsparseMatrix,nMatrix-method}
\alias{coerce,lsparseMatrix,nsparseMatrix-method}
\alias{coerce,matrix,lsparseMatrix-method}
\alias{coerce,numLike,lsparseMatrix-method}
\alias{which,lsparseMatrix-method}
% lgC
\alias{Arith,lgCMatrix,lgCMatrix-method}
\alias{Logic,lgCMatrix,lgCMatrix-method}
\alias{coerce,lgCMatrix,dgCMatrix-method}
\alias{coerce,lgCMatrix,lgRMatrix-method}
\alias{coerce,lgCMatrix,lgTMatrix-method}
\alias{coerce,lgCMatrix,lgeMatrix-method}
\alias{coerce,lgCMatrix,lsCMatrix-method}
\alias{coerce,lgCMatrix,lsTMatrix-method}
\alias{coerce,lgCMatrix,ltCMatrix-method}
\alias{coerce,lgCMatrix,ltTMatrix-method}
\alias{coerce,lgCMatrix,ngCMatrix-method}
\alias{coerce,matrix,lgCMatrix-method}
\alias{coerce,numLike,lgCMatrix-method}
% lgR
\alias{coerce,lgRMatrix,dgRMatrix-method}
\alias{coerce,lgRMatrix,lgCMatrix-method}
\alias{coerce,lgRMatrix,lgTMatrix-method}
\alias{coerce,lgRMatrix,lgeMatrix-method}
\alias{coerce,lgRMatrix,lsRMatrix-method}
\alias{coerce,lgRMatrix,ltRMatrix-method}
\alias{coerce,lgRMatrix,ngRMatrix-method}
\alias{coerce,matrix,lgRMatrix-method}
\alias{coerce,numLike,lgRMatrix-method}
% lgT
\alias{Arith,lgTMatrix,lgTMatrix-method}
\alias{Logic,lgTMatrix,lgTMatrix-method}
\alias{coerce,lgTMatrix,dgTMatrix-method}
\alias{coerce,lgTMatrix,lgCMatrix-method}
\alias{coerce,lgTMatrix,lgRMatrix-method}
\alias{coerce,lgTMatrix,lgeMatrix-method}
\alias{coerce,lgTMatrix,lsCMatrix-method}
\alias{coerce,lgTMatrix,lsTMatrix-method}
\alias{coerce,lgTMatrix,ltCMatrix-method}
\alias{coerce,lgTMatrix,ltTMatrix-method}
\alias{coerce,lgTMatrix,ngTMatrix-method}
\alias{coerce,matrix,lgTMatrix-method}
\alias{coerce,numLike,lgTMatrix-method}
\alias{which,lgTMatrix-method}
% ltC
\alias{Logic,ltCMatrix,ltCMatrix-method}
\alias{coerce,ltCMatrix,dtCMatrix-method}
\alias{coerce,ltCMatrix,lgCMatrix-method}
\alias{coerce,ltCMatrix,lgTMatrix-method}
\alias{coerce,ltCMatrix,lsCMatrix-method}
\alias{coerce,ltCMatrix,ltRMatrix-method}
\alias{coerce,ltCMatrix,ltTMatrix-method}
\alias{coerce,ltCMatrix,ltpMatrix-method}
\alias{coerce,ltCMatrix,ltrMatrix-method}
\alias{coerce,ltCMatrix,ntCMatrix-method}
\alias{coerce,matrix,ltCMatrix-method}
% ltR
\alias{coerce,ltRMatrix,dtRMatrix-method}
\alias{coerce,ltRMatrix,lgRMatrix-method}
\alias{coerce,ltRMatrix,lsRMatrix-method}
\alias{coerce,ltRMatrix,ltCMatrix-method}
\alias{coerce,ltRMatrix,ltTMatrix-method}
\alias{coerce,ltRMatrix,ltpMatrix-method}
\alias{coerce,ltRMatrix,ltrMatrix-method}
\alias{coerce,ltRMatrix,ntRMatrix-method}
\alias{coerce,matrix,ltRMatrix-method}
% ltT
\alias{coerce,ltTMatrix,dtTMatrix-method}
\alias{coerce,ltTMatrix,lgCMatrix-method}
\alias{coerce,ltTMatrix,lgTMatrix-method}
\alias{coerce,ltTMatrix,lsTMatrix-method}
\alias{coerce,ltTMatrix,ltCMatrix-method}
\alias{coerce,ltTMatrix,ltRMatrix-method}
\alias{coerce,ltTMatrix,ltpMatrix-method}
\alias{coerce,ltTMatrix,ltrMatrix-method}
\alias{coerce,ltTMatrix,ntTMatrix-method}
\alias{coerce,matrix,ltTMatrix-method}
\alias{which,ltTMatrix-method}
% lsC
\alias{Logic,lsCMatrix,lsCMatrix-method}
\alias{coerce,lsCMatrix,RsparseMatrix-method}
\alias{coerce,lsCMatrix,dsCMatrix-method}
\alias{coerce,lsCMatrix,lgCMatrix-method}
\alias{coerce,lsCMatrix,lgTMatrix-method}
\alias{coerce,lsCMatrix,lsRMatrix-method}
\alias{coerce,lsCMatrix,lsTMatrix-method}
\alias{coerce,lsCMatrix,lspMatrix-method}
\alias{coerce,lsCMatrix,lsyMatrix-method}
\alias{coerce,lsCMatrix,ltCMatrix-method}
\alias{coerce,lsCMatrix,nsCMatrix-method}
\alias{coerce,matrix,lsCMatrix-method}
% lsR
\alias{coerce,lsRMatrix,CsparseMatrix-method}
\alias{coerce,lsRMatrix,dsRMatrix-method}
\alias{coerce,lsRMatrix,lgRMatrix-method}
\alias{coerce,lsRMatrix,lsCMatrix-method}
\alias{coerce,lsRMatrix,lsTMatrix-method}
\alias{coerce,lsRMatrix,lspMatrix-method}
\alias{coerce,lsRMatrix,lsyMatrix-method}
\alias{coerce,lsRMatrix,ltRMatrix-method}
\alias{coerce,lsRMatrix,nsRMatrix-method}
\alias{coerce,matrix,lsRMatrix-method}
% lsT
\alias{coerce,lsTMatrix,dsTMatrix-method}
\alias{coerce,lsTMatrix,lgCMatrix-method}
\alias{coerce,lsTMatrix,lgTMatrix-method}
\alias{coerce,lsTMatrix,lsCMatrix-method}
\alias{coerce,lsTMatrix,lsRMatrix-method}
\alias{coerce,lsTMatrix,lspMatrix-method}
\alias{coerce,lsTMatrix,lsyMatrix-method}
\alias{coerce,lsTMatrix,ltTMatrix-method}
\alias{coerce,lsTMatrix,nsTMatrix-method}
\alias{coerce,matrix,lsTMatrix-method}
\alias{which,lsTMatrix-method}
%
\description{The \code{lsparseMatrix} class is a virtual class of sparse
  matrices with \code{TRUE}/\code{FALSE} or \code{NA} entries.  Only the
  positions of the elements that are \code{TRUE} are stored.

  These can be stored in the \dQuote{triplet} form (class
  \code{\linkS4class{TsparseMatrix}}, subclasses \code{lgTMatrix},
  \code{lsTMatrix}, and \code{ltTMatrix}) or in compressed
  column-oriented form (class \code{\linkS4class{CsparseMatrix}},
  subclasses \code{lgCMatrix}, \code{lsCMatrix}, and \code{ltCMatrix})
  or--\emph{rarely}--in compressed row-oriented form (class
  \code{\linkS4class{RsparseMatrix}}, subclasses \code{lgRMatrix},
  \code{lsRMatrix}, and \code{ltRMatrix}).  The second letter in the
  name of these non-virtual classes indicates \code{g}eneral,
  \code{s}ymmetric, or \code{t}riangular.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("lgCMatrix",
    ...)} and so on.  More frequently objects are created by coercion of
  a numeric sparse matrix to the logical form, e.g. in an expression
  \code{x != 0}.

  The logical form is also used in the symbolic analysis phase
  of an algorithm involving sparse matrices.  Such algorithms often
  involve two phases: a symbolic phase wherein the positions of the
  non-zeros in the result are determined and a numeric phase wherein the
  actual results are calculated.  During the symbolic phase only the
  positions of the non-zero elements in any operands are of interest,
  hence any numeric sparse matrices can be treated as logical sparse
  matrices.
}
\details{
  Note that triplet stored (\code{\linkS4class{TsparseMatrix}}) matrices
  such as \code{lgTMatrix} may contain duplicated pairs of indices
  \eqn{(i,j)} as for the corresponding numeric class
  \code{\linkS4class{dgTMatrix}} where for such pairs, the corresponding
  \code{x} slot entries are added.  For logical matrices, the \code{x}
  entries corresponding to duplicated index pairs \eqn{(i,j)} are
  \dQuote{added} as well if the addition is defined as logical \eqn{or},
  i.e., \dQuote{\code{TRUE + TRUE  |->  TRUE}} and
  \dQuote{\code{TRUE + FALSE  |->  TRUE}}.
  Note the use of \code{\link{uniqTsparse}()} for getting an internally
  unique representation without duplicated \eqn{(i,j)} entries.
}
\section{Slots}{
  \describe{
    \item{\code{x}:}{Object of class \code{"logical"}, i.e., either
      \code{TRUE}, \code{\link{NA}}, or \code{FALSE}.}
    \item{\code{uplo}:}{Object of class \code{"character"}. Must be
      either "U", for upper triangular, and "L", for lower
      triangular. Present in the triangular and symmetric classes but not
      in the general class.}
    \item{\code{diag}:}{Object of class \code{"character"}. Must be
      either \code{"U"}, for unit triangular (diagonal is all ones), or
      \code{"N"} for non-unit.  The implicit diagonal elements are not
      explicitly stored when \code{diag} is \code{"U"}.  Present in the
      triangular classes only.}
    \item{\code{p}:}{Object of class \code{"integer"} of pointers, one
      for each column (row), to the initial (zero-based) index of elements in
      the column.  Present in compressed column-oriented and compressed
      row-oriented forms only.}
    \item{\code{i}:}{Object of class \code{"integer"} of length nnzero
      (number of non-zero elements).  These are the row numbers for
      each TRUE element in the matrix.  All other elements are FALSE.
      Present in triplet and compressed column-oriented forms only.}
    \item{\code{j}:}{Object of class \code{"integer"} of length nnzero
      (number of non-zero elements).  These are the column numbers for
      each TRUE element in the matrix.  All other elements are FALSE.
      Present in triplet and compressed row-oriented forms only.}
    \item{\code{Dim}:}{Object of class \code{"integer"} - the dimensions
      of the matrix.}
  }
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "dgCMatrix", to = "lgCMatrix")}}
    \item{t}{\code{signature(x = "lgCMatrix")}: returns the transpose
      of \code{x}}

    \item{which}{\code{signature(x = "lsparseMatrix")}, semantically
      equivalent to \pkg{base} function \code{\link{which}(x, arr.ind)};
      for details, see the \code{\linkS4class{lMatrix}} class documentation.}
  }
}
\seealso{
  the class \code{\linkS4class{dgCMatrix}} and \code{\linkS4class{dgTMatrix}}
}
\examples{
(m <- Matrix(c(0,0,2:0), 3,5, dimnames=list(LETTERS[1:3],NULL)))
(lm <- (m > 1)) # lgC
!lm     # no longer sparse
stopifnot(is(lm,"lsparseMatrix"),
          identical(!lm, m <= 1))

data(KNex)
str(mmG.1 <- (KNex $ mm) > 0.1)# "lgC..."
table(mmG.1@x)# however with many ``non-structural zeros''
## from logical to nz_pattern -- okay when there are no NA's :
nmG.1 <- as(mmG.1, "nMatrix") # <<< has "TRUE" also where mmG.1 had FALSE
## from logical to "double"
dmG.1 <- as(mmG.1, "dMatrix") # has '0' and back:
lmG.1 <- as(dmG.1, "lMatrix")
stopifnot(identical(nmG.1, as((KNex $ mm) != 0,"nMatrix")),
          validObject(lmG.1),
          identical(lmG.1, mmG.1))

class(xnx <- crossprod(nmG.1))# "nsC.."
class(xlx <- crossprod(mmG.1))# "dsC.." : numeric
is0 <- (xlx == 0)
mean(as.vector(is0))# 99.3\% zeros: quite sparse, but
table(xlx@x == 0)# more than half of the entries are (non-structural!) 0
stopifnot(isSymmetric(xlx), isSymmetric(xnx),
          ## compare xnx and xlx : have the *same* non-structural 0s :
          sapply(slotNames(xnx),
                 function(n) identical(slot(xnx, n), slot(xlx, n))))
}
\keyword{classes}
\keyword{algebra}
