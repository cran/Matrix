\name{diagonalMatrix-class}
\title{Class "diagonalMatrix" of Diagonal Matrices}
\docType{class}
\alias{diagonalMatrix-class}
%
\alias{Math,diagonalMatrix-method}
\alias{Ops,diagonalMatrix,triangularMatrix-method}
\alias{as.logical,diagonalMatrix-method}
\alias{as.numeric,diagonalMatrix-method}
\alias{as.vector,diagonalMatrix-method}
\alias{cbind2,diagonalMatrix,sparseMatrix-method}
\alias{coerce,diagonalMatrix,CsparseMatrix-method}
\alias{coerce,diagonalMatrix,RsparseMatrix-method}
\alias{coerce,diagonalMatrix,TsparseMatrix-method}
\alias{coerce,diagonalMatrix,dMatrix-method}
\alias{coerce,diagonalMatrix,ddenseMatrix-method}
\alias{coerce,diagonalMatrix,denseMatrix-method}
\alias{coerce,diagonalMatrix,dsparseMatrix-method}
\alias{coerce,diagonalMatrix,generalMatrix-method}
\alias{coerce,diagonalMatrix,lMatrix-method}
\alias{coerce,diagonalMatrix,ldenseMatrix-method}
\alias{coerce,diagonalMatrix,lsparseMatrix-method}
\alias{coerce,diagonalMatrix,matrix-method}
\alias{coerce,diagonalMatrix,nMatrix-method}
\alias{coerce,diagonalMatrix,ndenseMatrix-method}
\alias{coerce,diagonalMatrix,nsparseMatrix-method}
\alias{coerce,diagonalMatrix,packedMatrix-method}
\alias{coerce,diagonalMatrix,sparseVector-method}
\alias{coerce,diagonalMatrix,nMatrix-method}
\alias{coerce,diagonalMatrix,symmetricMatrix-method}
\alias{coerce,diagonalMatrix,triangularMatrix-method}
\alias{coerce,diagonalMatrix,unpackedMatrix-method}
\alias{coerce,diagonalMatrix,vector-method}
\alias{coerce,matrix,diagonalMatrix-method}
\alias{determinant,diagonalMatrix,logical-method}
\alias{diag,diagonalMatrix-method}
\alias{diag<-,diagonalMatrix-method}
\alias{log,diagonalMatrix-method}
\alias{print,diagonalMatrix-method}
\alias{rbind2,diagonalMatrix,sparseMatrix-method}
\alias{show,diagonalMatrix-method}
\alias{summary,diagonalMatrix-method}
\alias{t,diagonalMatrix-method}
%
\description{
  Class "diagonalMatrix" is the virtual class of all diagonal matrices.
}
\section{Objects from the Class}{A virtual Class: No objects may be
  created from it.}
\section{Slots}{
  \describe{
    \item{\code{diag}:}{code{"character"} string, either \code{"U"} or
      \code{"N"}, where \code{"U"} means \sQuote{unit-diagonal}.}
    \item{\code{Dim}:}{matrix dimension, and}
    \item{\code{Dimnames}:}{the \code{\link{dimnames}}, a
      \code{\link{list}}, see the \code{\linkS4class{Matrix}} class
      description.  Typically \code{list(NULL,NULL)} for diagonal matrices.}
  }
}
\section{Extends}{
  Class \code{"\linkS4class{sparseMatrix}"}, directly.
}
\section{Methods}{
  These are just a subset of the signature for which defined methods.
  Currently, there are (too) many explicit methods defined in order to
  ensure efficient methods for diagonal matrices.

  \describe{
    \item{coerce}{\code{signature(from = "matrix", to = "diagonalMatrix")}: ... }
    \item{coerce}{\code{signature(from = "Matrix", to = "diagonalMatrix")}: ... }
    \item{coerce}{\code{signature(from = "diagonalMatrix", to = "generalMatrix")}: ... }
    \item{coerce}{\code{signature(from = "diagonalMatrix", to = "triangularMatrix")}: ... }
    \item{coerce}{\code{signature(from = "diagonalMatrix", to = "nMatrix")}: ... }
    \item{coerce}{\code{signature(from = "diagonalMatrix", to = "matrix")}: ... }
    \item{coerce}{\code{signature(from = "diagonalMatrix", to = "sparseVector")}: ... }
    \item{t}{\code{signature(x = "diagonalMatrix")}: ... }
    \cr and many more methods

    \item{solve}{\code{signature(a = "diagonalMatrix", b, ...)}: is
      trivially implemented, of course; see also \code{\link{solve-methods}}.}

    \item{which}{\code{signature(x = "nMatrix")}, semantically
      equivalent to \pkg{base} function \code{\link{which}(x, arr.ind)}.}

    \item{"Math"}{\code{signature(x = "diagonalMatrix")}: all these
      group methods return a \code{"diagonalMatrix"}, apart from
      \code{\link{cumsum}()} etc which return a \emph{vector} also for
      \pkg{base} \code{\link{matrix}}.}

    \item{*}{\code{signature(e1 = "ddiMatrix", e2="denseMatrix")}:
      arithmetic and other operators from the \code{\link[=S4groupGeneric]{Ops}}
      group have a few dozen explicit method definitions, in order to
      keep the results \emph{diagonal} in many cases, including the following:}
    \item{/}{\code{signature(e1 = "ddiMatrix", e2="denseMatrix")}:
      the result is from class \code{\linkS4class{ddiMatrix}} which is
      typically very desirable.  Note that when \code{e2} contains
      off-diagonal zeros or \code{\link{NA}}s, we implicitly use \eqn{0 / x = 0}, hence
      differing from traditional \R arithmetic (where \eqn{0 / 0
	\mapsto \mbox{NaN}}{0/0 |-> NaN}), in order to preserve sparsity.}

    \item{summary}{\code{(object = "diagonalMatrix")}: Returns
      an object of S3 class \code{"diagSummary"} which is the summary of
      the vector \code{object@x} plus a simple heading, and an
      appropriate \code{\link{print}} method.}
  }
}
\seealso{
  \code{\link{Diagonal}()} as constructor of these matrices, and
  \code{\link{isDiagonal}}.
  \code{\linkS4class{ddiMatrix}} and \code{\linkS4class{ldiMatrix}} are
  \dQuote{actual} classes extending \code{"diagonalMatrix"}.
}
\examples{
I5 <- Diagonal(5)
D5 <- Diagonal(x = 10*(1:5))
## trivial (but explicitly defined) methods:
stopifnot(identical(crossprod(I5), I5),
          identical(tcrossprod(I5), I5),
          identical(crossprod(I5, D5), D5),
          identical(tcrossprod(D5, I5), D5),
          identical(solve(D5), solve(D5, I5)),
          all.equal(D5, solve(solve(D5)), tolerance = 1e-12)
          )
solve(D5)# efficient as is diagonal

# an unusual way to construct a band matrix:
rbind2(cbind2(I5, D5),
       cbind2(D5, I5))
}
\keyword{classes}
