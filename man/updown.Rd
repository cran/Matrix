\name{updown-methods}
\title{Updating and Downdating Sparse Cholesky Factorizations}
%
\docType{methods}
\keyword{algebra}
\keyword{array}
\keyword{methods}
%
\alias{updown}
\alias{updown-methods}
%
\alias{updown,character,Matrix,CHMfactor-method}
\alias{updown,character,matrix,CHMfactor-method}
\alias{updown,logical,Matrix,CHMfactor-method}
\alias{updown,logical,dgCMatrix,CHMfactor-method}
\alias{updown,logical,dsCMatrix,CHMfactor-method}
\alias{updown,logical,dtCMatrix,CHMfactor-method}
\alias{updown,logical,matrix,CHMfactor-method}
%
\description{
  Computes a rank-\eqn{k} update or downdate of a sparse Cholesky
  factorization
  \deqn{P_{1} A P_{1}' = L_{1} D L_{1}' = L L'}{P1 * A * P1' = L1 * D * L1' = L * L'}
  which for some \eqn{k}-column matrix \eqn{C} is the factorization
  \deqn{P_{1} (A + s C C') P_{1}' = \tilde{L}_{1} \tilde{D} \tilde{L}_{1}' = \tilde{L} \tilde{L}'}{P1 * (A + s * C * C') * P1' = L~1 * D~ * L~1' = L~ * L~'}
  Here, \eqn{s = 1} for an update and \eqn{s = -1} for a downdate.
}
\usage{
updown(update, C, L)
}
\arguments{
  \item{update}{a logical (\code{TRUE} or \code{FALSE}) or
    character (\code{"+"} or \code{"-"}) indicating if \code{L}
    should be updated (or otherwise downdated).}
  \item{C}{a \link[=is.finite]{finite} matrix or
    \code{\linkS4class{Matrix}} such that
    \code{\link{tcrossprod}(C)} has the dimensions of \code{L}.}
  \item{L}{an object of class \code{\linkS4class{dCHMsimpl}} or
    \code{\linkS4class{dCHMsuper}} specifying a sparse Cholesky
    factorization.}
}
\value{
  A sparse Cholesky factorization with dimensions matching \code{L},
  typically of class \code{\linkS4class{dCHMsimpl}}.
}
\seealso{
  Classes
  \code{\linkS4class{dCHMsimpl}} and \code{\linkS4class{dCHMsuper}}
  and their methods, notably for generic function \code{\link{update}},
  which is \emph{not} equivalent to \code{updown(update = TRUE)}.

  Generic function \code{\link{Cholesky}}.
}
\author{Initial implementation by Nicholas Nagle, University of Tennessee.}
\references{
  Davis, T. A., Hager, W. W. (2001).
  Multiple-rank modifications of a sparse Cholesky factorization.
  \emph{SIAM Journal on Matrix Analysis and Applications},
  \emph{22}(4), 997-1013.
  \doi{10.1137/S0895479899357346}
}
\examples{
m <- sparseMatrix(i = c(3, 1, 3:2, 2:1), p = c(0:2, 4, 4, 6), x = 1:6,
                  dimnames = list(LETTERS[1:3], letters[1:5]))
uc0 <- Cholesky(A <- crossprod(m) + Diagonal(5))
uc1 <- updown("+", Diagonal(5, 1), uc0)
uc2 <- updown("-", Diagonal(5, 1), uc1)
stopifnot(all.equal(uc0, uc2))
\dontshow{
if(FALSE) {
## Hmm: this loses positive definiteness:
uc2 <- updown("-", Diagonal(5, 2), uc0)
image(show(as(uc0, "CsparseMatrix")))
image(show(as(uc2, "CsparseMatrix"))) # severely negative entries
}
}
}
