\name{isSymmetric-methods}
\title{Methods for Function 'isSymmetric' in Package 'Matrix'}
%
\docType{methods}
\keyword{array}
\keyword{programming}
\keyword{methods}
%
\alias{isSymmetric}
\alias{isSymmetric-methods}
%
\alias{isSymmetric,CsparseMatrix-method}
\alias{isSymmetric,RsparseMatrix-method}
\alias{isSymmetric,TsparseMatrix-method}
\alias{isSymmetric,denseMatrix-method}
\alias{isSymmetric,diagonalMatrix-method}
\alias{isSymmetric,indMatrix-method}
% tolerating numerical fuzz for [dz]Matrix:
\alias{isSymmetric,dgCMatrix-method}
\alias{isSymmetric,dgRMatrix-method}
\alias{isSymmetric,dgTMatrix-method}
\alias{isSymmetric,dgeMatrix-method}
\alias{isSymmetric,dtCMatrix-method}
\alias{isSymmetric,dtRMatrix-method}
\alias{isSymmetric,dtTMatrix-method}
\alias{isSymmetric,dtpMatrix-method}
\alias{isSymmetric,dtrMatrix-method}
%
\description{
  \code{isSymmetric} tests whether its argument is a symmetric square
  matrix, by default tolerating some numerical fuzz and requiring
  symmetric \code{[dD]imnames} in addition to symmetry in the
  mathematical sense.  \code{isSymmetric} is a generic function in
  \pkg{base}, which has a \link[=isSymmetric]{method} for traditional
  matrices of implicit \code{\link{class}} \code{"\link{matrix}"}.
  Methods are defined here for various proper and virtual classes
  in \pkg{Matrix}, so that \code{isSymmetric} works for all objects
  inheriting from virtual class \code{"\linkS4class{Matrix}"}.
}
\usage{
\S4method{isSymmetric}{denseMatrix}(object, checkDN = TRUE, \dots)
\S4method{isSymmetric}{CsparseMatrix}(object, checkDN = TRUE, \dots)
\S4method{isSymmetric}{RsparseMatrix}(object, checkDN = TRUE, \dots)
\S4method{isSymmetric}{TsparseMatrix}(object, checkDN = TRUE, \dots)
\S4method{isSymmetric}{diagonalMatrix}(object, checkDN = TRUE, \dots)
\S4method{isSymmetric}{indMatrix}(object, checkDN = TRUE, \dots)
\S4method{isSymmetric}{dgeMatrix}(object, checkDN = TRUE, tol = 100 * .Machine$double.eps, tol1 = 8 * tol, \dots)
\S4method{isSymmetric}{dgCMatrix}(object, checkDN = TRUE, tol = 100 * .Machine$double.eps, \dots)
}
\arguments{
  \item{object}{a \code{"Matrix"}.}
  \item{checkDN}{a \link{logical} indicating whether symmetry of the
    \code{Dimnames} \link{slot} of \code{object} should be checked.}
  \item{tol, tol1}{numerical tolerances allowing \emph{approximate}
    symmetry of numeric (rather than logical) matrices.  See also
    \code{\link{isSymmetric.matrix}}.}
  \item{\dots}{further arguments passed to methods
    (typically methods for \code{\link{all.equal}}).}
}
\value{
  A \link{logical}, either \code{TRUE} or \code{FALSE}
  (never \code{\link{NA}}).
}
\details{
  The \code{Dimnames} \link{slot} of \code{object}, say \code{dn},
  is considered to be symmetric if and only if
  \itemize{
    \item \code{dn[[1]]} and \code{dn[[2]]} are identical \emph{or}
    one is \code{NULL}; \emph{and}
    \item \code{ndn <- names(dn)} is \code{NULL} \emph{or}
    \code{ndn[1]} and \code{ndn[2]} are identical \emph{or}
    one is the empty string \code{""}.
  }
  Hence \code{list(a=nms, a=nms)} is considered to be \emph{symmetric},
  and so too are \code{list(a=nms, NULL)} and \code{list(NULL, a=nms)}.

  Note that this definition is \emph{looser} than that employed by
  \code{\link{isSymmetric.matrix}}, which requires \code{dn[1]} and
  \code{dn[2]} to be identical, where \code{dn} is the \code{dimnames}
  \link[=attr]{attribute} of a traditional matrix.
}
\seealso{
  \code{\link{forceSymmetric}};
  \code{\link{symmpart}} and \code{\link{skewpart}};
  virtual class \code{"\linkS4class{symmetricMatrix}"} and its subclasses.
}
\examples{
isSymmetric(Diagonal(4)) # TRUE of course
M <- Matrix(c(1,2,2,1), 2,2)
isSymmetric(M) # TRUE (*and* of formal class "dsyMatrix")
isSymmetric(as(M, "generalMatrix")) # still symmetric, even if not "formally"
isSymmetric(triu(M)) # FALSE

## Look at implementations:
showMethods("isSymmetric", includeDefs = TRUE) # includes S3 generic from base
}
