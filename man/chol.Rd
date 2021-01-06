\name{chol}
\title{Choleski Decomposition - 'Matrix' S4 Generic and Methods}
\alias{chol}
\alias{chol-methods}
\alias{chol,Matrix-method}
% dense methods
\alias{chol,ddenseMatrix-method}
\alias{chol,dgeMatrix-method}
\alias{chol,dpoMatrix-method}
\alias{chol,dppMatrix-method}
%   diagonal
\alias{chol,ddiMatrix-method}
\alias{chol,ldiMatrix-method}
% sparse methods
\alias{chol,dsparseMatrix-method}
\alias{chol,dsCMatrix-method}
\alias{chol,dsCMatrix-method}
\alias{chol,lsCMatrix-method}
\alias{chol,lsCMatrix-method}
\alias{chol,nsCMatrix-method}
\alias{chol,nsCMatrix-method}
\alias{chol,dtCMatrix-method}
\alias{chol,dtRMatrix-method}
\alias{chol,dtTMatrix-method}
%
\description{
  Compute the Choleski factorization of a real symmetric
  positive-definite square matrix.
}
\usage{
chol(x, \dots)
\S4method{chol}{dsCMatrix}(x, pivot = FALSE, \dots)
\S4method{chol}{dsparseMatrix}(x, pivot = FALSE, cache = TRUE, \dots)
}
\arguments{
  \item{x}{a (sparse or dense) square matrix, here inheriting from class
    \code{\linkS4class{Matrix}}; if \code{x} is not positive definite,
    an error is signalled.}
  \item{pivot}{logical indicating if pivoting is to be used.  Currently,
    this is \emph{not} made use of for dense matrices.}
  \item{cache}{logical indicating if the result should be cached in
    \code{x@factors}; note that this argument is experimental and only
    available for some sparse matrices.}
  \item{\dots}{potentially further arguments passed to methods.}
}
\details{
  Note that these Cholesky factorizations are typically \emph{cached} with
  \code{x} currently, and these caches are available in
  \code{x@factors}, which may be useful for the sparse case when
  \code{pivot = TRUE}, where the permutation can be retrieved; see also
  the examples.

  However, this should not be considered part of the API and made
  use of.  Rather consider \code{\link{Cholesky}()} in such situations,
  since \code{chol(x, pivot=TRUE)} uses the same algorithm (but not the
  same return value!) as \code{\link{Cholesky}(x, LDL=FALSE)} and
  \code{chol(x)} corresponds to
  \code{\link{Cholesky}(x, perm=FALSE, LDL=FALSE)}.
}
\section{Methods}{
  Use \code{\link{showMethods}(chol)} to see all; some are worth
  mentioning here:
  \describe{
    \item{chol}{\code{signature(x = "dgeMatrix")}: works via
      \code{"dpoMatrix"}, see class \code{\linkS4class{dpoMatrix}}.}

    \item{chol}{\code{signature(x = "dpoMatrix")}:
      Returns (and stores) the Cholesky decomposition of \code{x}, via
      LAPACK routines \code{dlacpy} and \code{dpotrf}.}

    \item{chol}{\code{signature(x = "dppMatrix")}:
      Returns (and stores) the Cholesky decomposition via LAPACK routine
      \code{dpptrf}.}

    \item{chol}{\code{signature(x = "dsCMatrix", pivot = "logical")}:
      Returns (and stores) the Cholesky decomposition of \code{x}.  If
      \code{pivot} is true, the Approximate Minimal Degree (AMD)
      algorithm is used to create a reordering of the rows and columns
      of \code{x} so as to reduce fill-in.}
  }
}
\value{
  a matrix of class \code{\linkS4class{Cholesky}},
  i.e., upper triangular: \eqn{R} such that \eqn{R'R = x} (if
      \code{pivot=FALSE}) \emph{or} \eqn{P' R'R P = x} (if
      \code{pivot=TRUE} and \eqn{P} is the corresponding permutation matrix).
}
\references{
  Timothy A. Davis (2006)
  \emph{Direct Methods for Sparse Linear Systems}, SIAM Series
  \dQuote{Fundamentals of Algorithms}.

  Tim Davis (1996),
  An approximate minimal degree ordering algorithm,
  \emph{SIAM J. Matrix Analysis and Applications},
  \bold{17}, 4, 886--905.
}
\seealso{The default from \pkg{base}, \code{\link[base]{chol}}; for more
  flexibility (but not returning a matrix!) \code{\link{Cholesky}}.
}
\examples{
showMethods(chol, inherited = FALSE) # show different methods

sy2 <- new("dsyMatrix", Dim = as.integer(c(2,2)), x = c(14, NA,32,77))
(c2 <- chol(sy2))#-> "Cholesky" matrix
stopifnot(all.equal(c2, chol(as(sy2, "dpoMatrix")), tolerance= 1e-13))
str(c2)

## An example where chol() can't work
(sy3 <- new("dsyMatrix", Dim = as.integer(c(2,2)), x = c(14, -1, 2, -7)))
try(chol(sy3)) # error, since it is not positive definite

## A sparse example --- exemplifying 'pivot'
(mm <- toeplitz(as(c(10, 0, 1, 0, 3), "sparseVector"))) # 5 x 5
(R <- chol(mm)) ## default:  pivot = FALSE
R2 <- chol(mm, pivot=FALSE)
stopifnot( identical(R, R2), all.equal(crossprod(R), mm) )
(R. <- chol(mm, pivot=TRUE))# nice band structure,
## but of course crossprod(R.) is *NOT* equal to mm
## --> see Cholesky() and its examples, for the pivot structure & factorization
stopifnot(all.equal(sqrt(det(mm)), det(R)),
          all.equal(prod(diag(R)), det(R)),
          all.equal(prod(diag(R.)), det(R)))

## a second, even sparser example:
(M2 <- toeplitz(as(c(1,.5, rep(0,12), -.1), "sparseVector")))
c2 <- chol(M2)
C2 <- chol(M2, pivot=TRUE)
## For the experts, check the caching of the factorizations:
ff <- M2@factors[["spdCholesky"]]
FF <- M2@factors[["sPdCholesky"]]
L1 <- as(ff, "Matrix")# pivot=FALSE: no perm.
L2 <- as(FF, "Matrix"); P2 <- as(FF, "pMatrix")
stopifnot(identical(t(L1), c2),
          all.equal(t(L2), C2, tolerance=0),#-- why not identical()?
          all.equal(M2, tcrossprod(L1)),             # M = LL'
          all.equal(M2, crossprod(crossprod(L2, P2)))# M = P'L L'P
         )
}
\keyword{algebra}
\keyword{array}
