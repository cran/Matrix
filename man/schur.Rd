\name{Schur}
\title{Schur Decomposition of a Matrix}
\usage{
Schur(x, vectors, \dots)
}
\alias{Schur}
\alias{Schur,dgeMatrix,logical-method}
\alias{Schur,dgeMatrix,missing-method}
\alias{Schur,ddenseMatrix,logical-method}
\alias{Schur,ddenseMatrix,missing-method}
\description{
    Computes the Schur decomposition and eigenvalues of a square matrix.
}
\arguments{
  \item{x}{
    numeric or complex square Matrix inheriting from class
    \code{"Matrix"}. Missing values (NAs) are not allowed.
  }
  \item{vectors}{logical.  When \code{TRUE} (the default), the Schur
    vectors are computed.
  }
  \item{\dots}{further arguments passed to or from other methods.}
}
\value{
  An object of class \code{c("schur.Matrix", "decomp")} whose
  attributes include the eigenvalues, the Schur quasi-triangular form
  of the matrix, and the Schur vectors (if requested).
}
\details{
  Based on the Lapack subroutine \code{dgees}.
}
\section{BACKGROUND}{
  If \code{A} is a square matrix, then \code{A = Q T t(Q)}, where
  \code{Q} is orthogonal, and \code{T} is upper quasi-triangular
  (nearly triangular with either 1 by 1 or 2 by 2 blocks on the
  diagonal).
  The eigenvalues of \code{A} are the same as those of \code{T},
  which are easy to compute.  The Schur form is used most often for
  computing non-symmetric eigenvalue decompositions, and for computing
  functions of matrices such as matrix exponentials.
}
\references{
  Anderson, E., et al. (1994).
  \emph{LAPACK User's Guide,}
  2nd edition, SIAM, Philadelphia.
}
\examples{
Schur(Hilbert(9))              # Schur factorization (real eigenvalues)

(A <- Matrix(round(rnorm(5*5, sd = 100)), nrow = 5))
(Sch.A <- Schur(A))
## The 'WR' and 'WI' components:
eTA <- eigen(Sch.A$T)
stopifnot(all.equal(eTA$values, eigen(A)$values, tol = 1e-13),
          all.equal(eTA$values,
                    local({z <- with(Sch.A, WR + 1i* WI)
                           z[order(Mod(z), decreasing=TRUE)]}), tol = 1e-13))
}
\keyword{algebra}
