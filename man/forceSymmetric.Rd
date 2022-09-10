\name{forceSymmetric}
\title{Force a Matrix to 'symmetricMatrix' Without Symmetry Checks}
%
\alias{forceSymmetric}
%
\alias{forceSymmetric,CsparseMatrix,character-method}
\alias{forceSymmetric,CsparseMatrix,missing-method}
\alias{forceSymmetric,RsparseMatrix,character-method}
\alias{forceSymmetric,RsparseMatrix,missing-method}
\alias{forceSymmetric,TsparseMatrix,character-method}
\alias{forceSymmetric,TsparseMatrix,missing-method}
\alias{forceSymmetric,diagonalMatrix,character-method}
\alias{forceSymmetric,diagonalMatrix,missing-method}
\alias{forceSymmetric,indMatrix,character-method}
\alias{forceSymmetric,indMatrix,missing-method}
\alias{forceSymmetric,matrix,character-method}
\alias{forceSymmetric,matrix,missing-method}
\alias{forceSymmetric,packedMatrix,character-method}
\alias{forceSymmetric,packedMatrix,missing-method}
\alias{forceSymmetric,unpackedMatrix,character-method}
\alias{forceSymmetric,unpackedMatrix,missing-method}
%
\description{
  Force a square matrix \code{x} to a \code{\linkS4class{symmetricMatrix}},
  \bold{without} a symmetry check as it would be applied for \code{as(x,
    "symmetricMatrix")}.
}
\usage{
forceSymmetric(x, uplo)
}
\arguments{
  \item{x}{any square matrix (of numbers), either \dQuote{"traditional"}
    (\code{\link{matrix}}) or inheriting from
    \code{\linkS4class{Matrix}}.}
  \item{uplo}{optional string, \code{"U"} or \code{"L"} indicating which
    \dQuote{triangle} half of \code{x} should determine the result.  The
    default is \code{"U"} unless \code{x} already has a \code{uplo} slot
    (i.e., when it is \code{\linkS4class{symmetricMatrix}}, or
    \code{\linkS4class{triangularMatrix}}), where the default will be
    \code{x@uplo}.}
}
% \details{
%
% }
\value{
  a square matrix inheriting from class
  \code{\linkS4class{symmetricMatrix}}.
}
\seealso{\code{\link{symmpart}} for the symmetric part of a matrix, or
  the coercions \code{as(x, <symmetricMatrix class>)}.
}
\examples{
 ## Hilbert matrix
 i <- 1:6
 h6 <- 1/outer(i - 1L, i, "+")
 sd <- sqrt(diag(h6))
 hh <- t(h6/sd)/sd # theoretically symmetric
 isSymmetric(hh, tol=0) # FALSE; hence
 try( as(hh, "symmetricMatrix") ) # fails, but this works fine:
 H6 <- forceSymmetric(hh)

 ## result can be pretty surprising:
 (M <- Matrix(1:36, 6))
 forceSymmetric(M) # symmetric, hence very different in lower triangle
 (tm <- tril(M))
 forceSymmetric(tm)
}
\keyword{array}
