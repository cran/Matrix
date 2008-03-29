\name{sparseLU-class}
\docType{class}
\alias{sparseLU-class}
\title{Sparse LU decomposition of a square sparse matrix}
\description{Objects of this class contain the components of the LU
  decomposition of a sparse square matrix.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("sparseLU",
    ...)} but are more commonly created by function \code{\link{lu}}
  applied to a sparse matrix, such as a matrix of class
  \code{\linkS4class{dgCMatrix}}.
}
\section{Slots}{
  \describe{
    \item{\code{L}:}{Object of class \code{"\linkS4class{dgCMatrix}"} The lower
      triangular factor from the left.}
    \item{\code{U}:}{Object of class \code{"\linkS4class{dgCMatrix}"} The upper
      triangular factor from the right.}
    \item{\code{p}:}{Object of class \code{"integer"} Permutation
      applied from the left. }
    \item{\code{q}:}{Object of class \code{"integer"} Permutation
      applied from the right.}
    \item{\code{Dim}:}{the dimension of the original matrix; inherited
      from class \code{\linkS4class{MatrixFactorization}}.}
  }
}
\section{Extends}{
  Class \code{"\linkS4class{LU}"}, directly.
  Class \code{"\linkS4class{MatrixFactorization}"}, by class \code{"LU"}.
}
\section{Methods}{
  No methods defined with class "sparseLU" in the signature.
}
%\references{}
%\author{}
\note{
  The decomposition is of the form
  \deqn{A = PLUQ}{A = PLUQ}
  where all matrices are sparse and of size \eqn{n\times n}{n by n}.
  The matrices \eqn{P} and \eqn{Q} are permutation matrices, \eqn{L} is
  lower triangular and \eqn{U} is upper triangular.
}
\seealso{
  \code{\link{lu}}, \code{\link[base]{solve}}, \code{\linkS4class{dgCMatrix}}
}
\examples{
  ## see  examples(lu)
}
\keyword{classes}
