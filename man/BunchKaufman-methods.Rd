\name{BunchKaufman-methods}
\title{Bunch-Kaufman Decomposition Methods}
\docType{methods}
\alias{BunchKaufman}
\alias{BunchKaufman-methods}
\alias{BunchKaufman,dspMatrix-method}
\alias{BunchKaufman,dsyMatrix-method}
\alias{BunchKaufman,matrix-method}
\description{
  The Bunch-Kaufman Decomposition of a square symmetric matrix \eqn{A}
  is \eqn{A = P LDL' P'} where \eqn{P} is a permutation matrix, \eqn{L}
  is \emph{unit}-lower triangular and \eqn{D} is \emph{block}-diagonal
  with blocks of dimension \eqn{1\times 1}{1 x 1} or \eqn{2\times2}{2 x 2}.

  This is generalization of a pivoting \eqn{LDL'} Cholesky decomposition.
}
\usage{
\S4method{BunchKaufman}{dsyMatrix}(x, \dots)
\S4method{BunchKaufman}{dspMatrix}(x, \dots)
\S4method{BunchKaufman}{matrix}(x, uplo = NULL, \dots)
}
\arguments{
  \item{x}{a symmetric square matrix.}
  \item{uplo}{optional string, \code{"U"} or \code{"L"} indicating which
    \dQuote{triangle} half of \code{x} should determine the result.  The
    default is \code{"U"} unless \code{x} has a \code{uplo} slot which is
    the case for those inheriting from class \code{\linkS4class{symmetricMatrix}},
    where \code{x@uplo} will be used.}
  \item{\dots}{potentially further arguments passed to methods.}
}
\section{Methods}{
  Currently, only methods for \bold{dense} numeric symmetric matrices
  are implemented.  To compute the Bunch-Kaufman decomposition, the methods
  use either one of two Lapack routines:
  \describe{
    \item{\code{x = "dspMatrix"}}{routine \code{dsptrf()}; whereas}
    \item{\code{x = "dsyMatrix"}}{, and}
    \item{\code{x = "matrix"}}{use \code{dsytrf()}.}
  }
}
\value{
  an object of class \code{\linkS4class{BunchKaufman}}, which can also
  be used as a (triangular) matrix directly.  Somewhat amazingly,
  it inherits its \code{uplo} slot from \code{x}.% FIXME?
}
\details{
  FIXME: We really need an \code{expand()} method in order to \emph{work} with the result!
}
\references{
  The original LAPACK source code, including documentation;
  \url{https://www.netlib.org/lapack/double/dsytrf.f}
  and
  \url{https://www.netlib.org/lapack/double/dsptrf.f}
}
\seealso{
  The resulting class, \code{\linkS4class{BunchKaufman}}.
  Related decompositions are the LU,  \code{\link{lu}}, and the Cholesky,
  \code{\link{chol}} (and for \emph{sparse} matrices,
  \code{\link{Cholesky}}).
}
\examples{
data(CAex)
dim(CAex)
isSymmetric(CAex)# TRUE
CAs <- as(CAex, "symmetricMatrix")
if(FALSE) # no method defined yet for *sparse* :
   bk. <- BunchKaufman(CAs)
## does apply to *dense* symmetric matrices:
bkCA <- BunchKaufman(as(CAs, "denseMatrix"))
bkCA

image(bkCA)# shows how sparse it is, too
str(R.CA <- as(bkCA, "sparseMatrix"))
## an upper triangular 72x72 matrix with only 144 non-zero entries
}
\keyword{methods}

