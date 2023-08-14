\name{drop0}
\title{Drop Non-Structural Zeros from a Sparse Matrix}
%
\keyword{array}
\keyword{manip}
\keyword{utilities}
%
\alias{drop0}
%
\description{
  Deletes \dQuote{non-structural} zeros (i.e., zeros stored explicitly,
  in memory) from a sparse matrix and returns the result.
}
\usage{
drop0(x, tol = 0, is.Csparse = NA, give.Csparse = TRUE)
}
\arguments{
  \item{x}{a \code{\linkS4class{Matrix}}, typically inheriting
    from virtual class \code{\linkS4class{sparseMatrix}}.
    \code{\linkS4class{denseMatrix}} and traditional vectors and
    matrices are coerced to \code{\linkS4class{CsparseMatrix}},
    with zeros dropped automatically, hence users passing such
    \code{x} should consider \code{as(x, "CsparseMatrix")} instead,
    notably in the \code{tol = 0} case.}
  \item{tol}{a non-negative number.  If \code{x} is numeric,
    then entries less than or equal to \code{tol} in absolute value
    are deleted.}
  \item{is.Csparse}{a logical used only if \code{give.Csparse}
    is \code{TRUE}, indicating if \code{x} already inherits from
    virtual class \code{\linkS4class{CsparseMatrix}}, in which
    case coercion is not attempted, permitting some (typically
    small) speed-up.}
  \item{give.Csparse}{a logical indicating if the result must
    inherit from virtual class \code{\linkS4class{CsparseMatrix}}.
    If \code{FALSE} and \code{x} inherits from
    \code{\linkS4class{RsparseMatrix}},
    \code{\linkS4class{TsparseMatrix}}, or
    \code{\linkS4class{indMatrix}},
    then the result preserves the class of \code{x}.
    The default value is \code{TRUE} only for backwards compatibility.}
}
\value{
  A \code{\linkS4class{sparseMatrix}}, the result of deleting
  non-structural zeros from \code{x}, possibly after coercion.
}
\note{
  \code{drop0} is sometimes called in conjunction with
  \code{\link{zapsmall}}, e.g., when dealing with sparse
  matrix products; see the example.
}
\seealso{
  Function \code{\link{sparseMatrix}}, for constructing objects
  inheriting from virtual class \code{\linkS4class{sparseMatrix}};
  \code{\link{nnzero}}.
}
\examples{
(m <- sparseMatrix(i = 1:8, j = 2:9, x = c(0:2, 3:-1),
                   dims = c(10L, 20L)))
drop0(m)

## A larger example:
t5 <- new("dtCMatrix", Dim = c(5L, 5L), uplo = "L",
          x = c(10, 1, 3, 10, 1, 10, 1, 10, 10),
          i = c(0L,2L,4L, 1L, 3L,2L,4L, 3L, 4L),
          p = c(0L, 3L, 5L, 7:9))
TT <- kronecker(t5, kronecker(kronecker(t5, t5), t5))
IT <- solve(TT)
I. <- TT \%*\% IT ;  nnzero(I.) # 697 ( == 625 + 72 )
I.0 <- drop0(zapsmall(I.))
## which actually can be more efficiently achieved by
I.. <- drop0(I., tol = 1e-15)
stopifnot(all(I.0 == Diagonal(625)), nnzero(I..) == 625)
}
