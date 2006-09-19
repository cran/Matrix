\name{pMatrix-class}
\docType{class}
\alias{pMatrix-class}
\alias{-,pMatrix,missing-method}
\alias{\%*\%,matrix,pMatrix-method}
\alias{\%*\%,pMatrix,matrix-method}
\alias{\%*\%,pMatrix,pMatrix-method}
\alias{\%*\%,Matrix,pMatrix-method}
\alias{\%*\%,pMatrix,Matrix-method}
\alias{coerce,integer,pMatrix-method}
\alias{coerce,numeric,pMatrix-method}
\alias{coerce,pMatrix,matrix-method}
\alias{coerce,pMatrix,ngTMatrix-method}
\alias{coerce,pMatrix,TsparseMatrix-method}
\alias{solve,pMatrix,missing-method}
\alias{t,pMatrix-method}
\title{Permutation matrices}
\description{The \code{"pMatrix"} class is the class of permutation
  matrices, stored as 1-based integer permutation vectors.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("pMatrix", ...)}
  or by coercion from an integer permutation vector, see below.
}
\section{Slots}{
  \describe{
    \item{\code{perm}:}{An integer, 1-based permutation vector, i.e.
      an integer vector of length \code{Dim[1]} whose elements form a
      permutation of \code{1:Dim[1]}.}
    \item{\code{Dim}:}{Object of class \code{"integer"}. The dimensions
      of the matrix which must be a two-element vector of equal,
      non-negative integers.}
    \item{\code{Dimnames}:}{list of length two; each component
      containing NULL or a \code{\link{character}} vector length
      equal the corresponding \code{Dim} element.}
  }
}
\section{Extends}{
  Class \code{"Matrix"}, directly.
}
\section{Methods}{
  \describe{
    \item{\%*\%}{\code{signature(x = "matrix", y = "pMatrix")} and other
      signatures (use \code{showMethods("\%*\%", class="pMatrix")}): ... }
    \item{coerce}{\code{signature(from = "integer", to = "pMatrix")}:
      This is enables typical \code{"pMatrix"} construction, given
      a permutation vector of \code{1:n}, see the first example.}
    \item{coerce}{\code{signature(from = "numeric", to = "pMatrix")}:
      a user convenience, to allow \code{as(perm, "pMatrix")} for
      numeric \code{perm} with integer values.}
    \item{coerce}{\code{signature(from = "pMatrix", to = "matrix")}: ... }
    \item{coerce}{\code{signature(from = "pMatrix", to = "ngTMatrix")}:
      coercion to sparse logical matrix of class \code{\linkS4class{ngTMatrix}}.}
    \item{solve}{\code{signature(a = "pMatrix", b = "missing")}: return
      the inverse permutation matrix.}
    \item{t}{\code{signature(x = "pMatrix")}: return the transpose of
      the permuation matrix (which is also the inverse of the
      permutation matrix).}
  }
}
\note{
  The inverse of the typical \code{"pMatrix"} constructor,
  \code{P <- as(ip, "pMatrix")} is simply \code{ip <- P@perm}.
}
%\seealso{}
\examples{
(pm1 <- as(as.integer(c(2,3,1)), "pMatrix"))
t(pm1) # is the same as
solve(pm1)
pm1 \%*\% t(pm1) # check that the transpose is the inverse
stopifnot(identical(diag(3), as(pm1 \%*\% t(pm1), "matrix")))

set.seed(11)
## random permutation matrix :
(p10 <- as(sample(10),"pMatrix"))

## Permute rows / columns of a numeric matrix :
(mm <- round(array(rnorm(3 * 3), c(3, 3)), 2))
mm \%*\% pm1
pm1 \%*\% mm
try(as(as.integer(c(3,3,1)), "pMatrix"))# Error: not a permutation

as(pm1, "ngTMatrix")
p10[1:7, 1:4] # gives an "ngTMatrix" (most economic!)
}
\keyword{classes}
