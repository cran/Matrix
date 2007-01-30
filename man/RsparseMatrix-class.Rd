\name{RsparseMatrix-class}
\docType{class}
\alias{RsparseMatrix-class}
%
\alias{coerce,RsparseMatrix,CsparseMatrix-method}
\alias{coerce,RsparseMatrix,TsparseMatrix-method}
\alias{coerce,RsparseMatrix,dgeMatrix-method}
\alias{coerce,RsparseMatrix,matrix-method}
\alias{coerce,RsparseMatrix,dsparseMatrix-method}
\alias{coerce,RsparseMatrix,lsparseMatrix-method}
\alias{coerce,RsparseMatrix,nsparseMatrix-method}
\alias{coerce,RsparseMatrix,dMatrix-method}
\alias{coerce,RsparseMatrix,lMatrix-method}
\alias{coerce,RsparseMatrix,nMatrix-method}
% all other ..RMatrix coercions -- here as well for now --
\alias{coerce,matrix,dgRMatrix-method}
%
\alias{coerce,dgCMatrix,dgRMatrix-method}
\alias{coerce,dgRMatrix,CsparseMatrix-method}
\alias{coerce,dgRMatrix,dgCMatrix-method}
\alias{coerce,dgRMatrix,dgTMatrix-method}
\alias{coerce,dgRMatrix,dgeMatrix-method}
%\alias{coerce,dgRMatrix,matrix-method}
\alias{coerce,dgTMatrix,dgRMatrix-method}
\alias{coerce,dgeMatrix,dgRMatrix-method}
\alias{coerce,dsRMatrix,dsCMatrix-method}
\alias{coerce,dsRMatrix,dsTMatrix-method}
\alias{coerce,dsRMatrix,dsyMatrix-method}
\alias{coerce,dtRMatrix,dtCMatrix-method}
\alias{coerce,dtRMatrix,dtTMatrix-method}
\alias{coerce,dtRMatrix,dtrMatrix-method}
\alias{coerce,lgRMatrix,lgCMatrix-method}
\alias{coerce,lgRMatrix,lgTMatrix-method}
\alias{coerce,lgRMatrix,lgeMatrix-method}
%\alias{coerce,lgRMatrix,matrix-method}
\alias{coerce,lsRMatrix,lsCMatrix-method}
\alias{coerce,lsRMatrix,lsTMatrix-method}
\alias{coerce,lsRMatrix,lsyMatrix-method}
\alias{coerce,ltRMatrix,ltCMatrix-method}
\alias{coerce,ltRMatrix,ltTMatrix-method}
\alias{coerce,ltRMatrix,ltrMatrix-method}
%\alias{coerce,ngRMatrix,matrix-method}
\alias{coerce,ngRMatrix,ngCMatrix-method}
\alias{coerce,ngRMatrix,ngTMatrix-method}
\alias{coerce,ngRMatrix,ngeMatrix-method}
\alias{coerce,nsRMatrix,nsCMatrix-method}
\alias{coerce,nsRMatrix,nsTMatrix-method}
\alias{coerce,nsRMatrix,nsyMatrix-method}
\alias{coerce,ntRMatrix,ntCMatrix-method}
\alias{coerce,ntRMatrix,ntTMatrix-method}
\alias{coerce,ntRMatrix,ntrMatrix-method}
%
\alias{t,RsparseMatrix-method}
%
\title{Class "RsparseMatrix" of Sparse Matrices in Column-compressed Form}
\description{The \code{"RsparseMatrix"} class is the virtual class of
  all sparse matrices coded in sorted compressed row-oriented form.
  Since it is a virtual class, no objects may be created from it.  See
  \code{showClass("RsparseMatrix")} for its subclasses.
}
\section{Slots}{
  \describe{
    \item{\code{j}:}{Object of class \code{"integer"} of length \code{nnzero}
      (number of non-zero elements).  These are the row numbers for
      each non-zero element in the matrix.}
    \item{\code{p}:}{Object of class \code{"integer"} of pointers, one
      for each row, to the initial (zero-based) index of elements in
      the row.}
    \item{\code{Dim}, \code{Dimnames}:}{inherited from
      the superclass, see \code{\linkS4class{sparseMatrix}}.}
  }
}
\section{Extends}{
  Class \code{"sparseMatrix"}, directly.
  Class \code{"Matrix"}, by class \code{"sparseMatrix"}.
}
\section{Methods}{
  Only \bold{few} methods are defined currently on purpose, since we
  rather use the \code{\linkS4class{CsparseMatrix}} in \pkg{Matrix}.
  \describe{
    \item{t}{\code{signature(x = "RsparseMatrix")}: ... }
    \item{coerce}{\code{signature(from = "RsparseMatrix", to = "CsparseMatrix")}: ... }
    \item{coerce}{\code{signature(from = "RsparseMatrix", to = "TsparseMatrix")}: ... }
  }
}
\seealso{
  its superclass, \code{\linkS4class{sparseMatrix}}, and, e.g., class
  \code{\linkS4class{dgRMatrix}} for the links to other classes.
}
\examples{
showClass("RsparseMatrix")
}
\keyword{classes}
