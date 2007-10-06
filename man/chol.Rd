\name{chol}
\title{Choleski Decomposition - 'Matrix' S4 Generic}
\alias{chol}
\alias{chol,Matrix,ANY-method}
% dense methods
\alias{chol,ddenseMatrix,ANY-method}
\alias{chol,dgeMatrix,ANY-method}
\alias{chol,dpoMatrix-method}
\alias{chol,dppMatrix-method}
\alias{chol,dpoMatrix,ANY-method}
\alias{chol,dppMatrix,ANY-method}
%   diagonal
\alias{chol,ddiMatrix,ANY-method}
\alias{chol,ldiMatrix,ANY-method}
% sparse methods
\alias{chol,dsparseMatrix,ANY-method}
\alias{chol,dsCMatrix,missing-method}
\alias{chol,dsCMatrix,logical-method}
\alias{chol,lsCMatrix,logical-method}
\alias{chol,lsCMatrix,missing-method}
\alias{chol,nsCMatrix,logical-method}
\alias{chol,nsCMatrix,missing-method}
%
\description{
  Compute the Choleski factorization of a real symmetric
  positive-definite square matrix.
}
\usage{
chol(x, pivot = FALSE, \dots)
}
\arguments{
  \item{x}{a (sparse or dense) square matrix, here inheriting from class
    \code{\linkS4class{Matrix}}; if \code{x} is not positive definite,
    an error is signalled.}
  \item{pivot}{logical indicating if pivoting is used.}
  \item{\dots}{potentially further arguments passed to methods.}
}
% \details{
% }
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
      \code{pivot} is \code{TRUE} (the default) Metis is used to create
      a reordering of the rows and columns of \code{x} so as to minimize
      fill-in.}

  }
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\seealso{The default from \pkg{base}, \code{\link[base]{chol}}.
}
\examples{
showMethods(chol, inherited = FALSE) # show different methods
}
\keyword{algebra}
\keyword{array}
