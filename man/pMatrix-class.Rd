\name{pMatrix-class}
\title{Permutation matrices}
%
\docType{class}
\keyword{array}
\keyword{classes}
%
\alias{pMatrix-class}
%
\alias{coerce,matrix,pMatrix-method}
\alias{coerce,numeric,pMatrix-method}
\alias{determinant,pMatrix,logical-method}
\alias{t,pMatrix-method}
%
\description{
  The \code{pMatrix} class is the class of \emph{permutation} matrices,
  stored as 1-based integer permutation vectors.  A permutation
  matrix is a square matrix whose rows \emph{and} columns are all
  standard unit vectors.  It follows that permutation matrices are
  a special case of \emph{index} matrices (hence \code{pMatrix}
  is defined as a direct subclass of \code{\linkS4class{indMatrix}}).
  
  Multiplying a matrix on the left by a permutation matrix is
  equivalent to permuting its rows.  Analogously, multiplying a
  matrix on the right by a permutation matrix is equivalent to
  permuting its columns.  Indeed, such products are implemented in
  \pkg{Matrix} as indexing operations; see \sQuote{Details} below.
}
\section{Objects from the Class}{
  Objects can be created explicitly with calls of the form
  \code{new("pMatrix", ...)}, but they are more commonly created
  by coercing 1-based integer index vectors, with calls of the
  form \code{as(., "pMatrix")}; see \sQuote{Methods} below.
}
\section{Slots}{
  \describe{
    \item{\code{margin},\code{perm}}{inherited from superclass
      \code{\linkS4class{indMatrix}}.  Here, \code{perm} is an
      integer vector of length \code{Dim[1]} and a permutation
      of \code{1:Dim[1]}}.
    \item{\code{Dim},\code{Dimnames}}{inherited from virtual
      superclass \code{\linkS4class{Matrix}}.}
  }
}
\section{Extends}{
  Class \code{"\linkS4class{indMatrix}"}, directly.
}
\section{Methods}{
  \describe{
    \item{\code{\%*\%}}{\code{signature(x = "pMatrix", y = "Matrix")}
      and others listed by \code{showMethods("\%*\%", classes = "pMatrix")}:
      matrix products implemented where appropriate as indexing operations.}
    \item{\code{coerce}}{\code{signature(from = "numeric", to = "pMatrix")}:
      supporting typical \code{pMatrix} construction from a vector
      of positive integers, specifically a permutation of \code{1:n}.
      Row permutation is assumed.}
    \item{t}{\code{signature(x = "pMatrix")}:
      the transpose, which is a \code{pMatrix} with identical
      \code{perm} but opposite \code{margin}.  Coincides with
      the inverse, as permutation matrices are orthogonal.}
    \item{solve}{\code{signature(a = "pMatrix", b = "missing")}:
      the inverse permutation matrix, which is a \code{pMatrix}
      with identical \code{perm} but opposite \code{margin}.
      Coincides with the transpose, as permutation matrices are
      orthogonal.  See \code{showMethods("solve", classes = "pMatrix")}
      for more signatures.}
    \item{determinant}{\code{signature(x = "pMatrix", logarithm = "logical")}:
      always returning 1 or -1, as permutation matrices are orthogonal.  
      In fact, the result is exactly the \emph{sign} of the permutation.}
  }
}
\details{
  By definition, a permutation matrix is both a row index matrix
  and a column index matrix.  However, the \code{perm} slot of
  a \code{pMatrix} cannot be used interchangeably as a row index
  vector and column index vector.  If \code{margin=1}, then
  \code{perm} is a row index vector, and the corresponding column
  index vector can be computed as \code{\link{invPerm}(perm)}, i.e.,
  by inverting the permutation.  Analogously, if \code{margin=2},
  then \code{perm} and \code{invPerm(perm)} are column and row
  index vectors, respectively.

  Given an \code{n}-by-\code{n} row permutation matrix \code{P}
  with \code{perm} slot \code{p} and a matrix \code{M} with
  conformable dimensions, we have
  \tabular{lclcl}{
    \eqn{P M} \tab = \tab \code{P \%*\% M}        \tab = \tab \code{M[p, ]}\cr
    \eqn{M P} \tab = \tab \code{M \%*\% P}        \tab = \tab \code{M[, i(p)]}\cr
    \eqn{P'M} \tab = \tab \code{crossprod(P, M)}  \tab = \tab \code{M[i(p), ]}\cr
    \eqn{MP'} \tab = \tab \code{tcrossprod(M, P)} \tab = \tab \code{M[, p]}\cr
    \eqn{P'P} \tab = \tab \code{crossprod(P)}     \tab = \tab \code{Diagonal(n)}\cr
    \eqn{PP'} \tab = \tab \code{tcrossprod(P)}    \tab = \tab \code{Diagonal(n)}
  }
  where \code{i := invPerm}.
}
\seealso{
  Superclass \code{\linkS4class{indMatrix}} of index matrices,
  for many inherited methods; \code{\link{invPerm}}, for computing
  inverse permutation vectors.
}
\examples{
\dontshow{ % for R_DEFAULT_PACKAGES=NULL
library(stats, pos = "package:base", verbose = FALSE)
}
(pm1 <- as(as.integer(c(2,3,1)), "pMatrix"))
t(pm1) # is the same as
solve(pm1)
pm1 \%*\% t(pm1) # check that the transpose is the inverse
stopifnot(all(diag(3) == as(pm1 \%*\% t(pm1), "matrix")),
          is.logical(as(pm1, "matrix")))

set.seed(11)
## random permutation matrix :
(p10 <- as(sample(10),"pMatrix"))

## Permute rows / columns of a numeric matrix :
(mm <- round(array(rnorm(3 * 3), c(3, 3)), 2))
mm \%*\% pm1
pm1 \%*\% mm
try(as(as.integer(c(3,3,1)), "pMatrix"))# Error: not a permutation

as(pm1, "TsparseMatrix")
p10[1:7, 1:4] # gives an "ngTMatrix" (most economic!)

## row-indexing of a <pMatrix> keeps it as an <indMatrix>:
p10[1:3, ]
}
