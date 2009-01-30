\name{spqr}
\title{SuiteSparse QR factorization of Sparse Matrix}
\alias{spqr}
\alias{spqr,CsparseMatrix-method}
\alias{spqr,TsparseMatrix-method}
\alias{spqr,dgCMatrix-method}
%
\description{
 SuiteSparse QR ("SPQR") implements a multifrontal sparse QR
 factorization method.
}
\usage{
spqr(x, ...)

\S4method{spqr}{dgCMatrix}(x, econ = 0,
     ordering = c("default", "fixed", "natural", "amd", "colamd",
                  "metis", "best", "bestamd"),
     tol = -2)
}
\arguments{
  \item{x}{a sparse matrix, i.e., an \R object inheriting from class
    \code{\linkS4class{sparseMatrix}}.}
  \item{econ}{integer for \dQuote{economy}, i.e. specifying the number
    of rows of \eqn{R} or columns of \eqn{Q} to return.  The default is
    \code{m <- nrow(x)}.}
  \item{ordering}{a string specifying the ordering strategy.}
  \item{tol}{a number specifying the tolerance to be applied. If
     \code{tol <= -2}, a default tolerance is used; if \eqn{-2 < tol <
     0}, then no tolerance is used, and if \eqn{tol >= 0} it is used as
   tolerance, namely columns with 2-norm <= tol are treated as 0.}
  \item{\dots}{further arguments passed to methods.}
}
% \details{
%   ~~ If necessary, more details than the description above ~~
% }
\value{
  an object of class \code{"\linkS4class{SPQR}"}, notably with slots
  \code{Q} and \code{R} containing the factorizing matrices.
}
\section{Methods}{
  \describe{
    \item{x = "dgCMatrix"}{ ~~describe this method here }
    \item{x = "CsparseMatrix"}{ ~~describe this method here }
    \item{x = "TsparseMatrix"}{ ~~describe this method here }
  }
}
\references{\url{http://www.cise.ufl.edu/research/sparse/SPQR/} has
  links to publications and a user manual.
}
\seealso{
  \code{\link[Matrix]{qr}} for \emph{long} matrices, uses a simpler
  algorithm from library CSparse.
}
\examples{
data(KNex); mm <- KNex$mm
if(.Platform$OS.type != "windows") # C code seg.faults on some Windows versions
str(mm.QR <- spqr(mm))
}
\keyword{algebra}
\keyword{array}
