\name{schur}
\alias{schur}
\title{
Schur Decomposition of a Matrix --- Generic Function
}
\description{
Computes the Schur decomposition (including eigenvalues) of a square matrix.
In S-PLUS, only one class has a method for this function -- the
\code{"Matrix"} class from the \code{Matrix} library. It is defined as a method to
allow users to easily incorporate customized versions into S-PLUS.
}
\usage{
schur(x, \dots)
}
\arguments{
\item{x}{
a square matrix. No missing values or IEEE special values are allowed.
}
\item{...}{
most methods will have additional arguments, for example an argument
indicating whether or not to compute the Schur vectors.
}}
\value{
a representation of all or part of the Schur decomposition of \code{x}.
}
\section{BACKGROUND}{
If \code{A} is a square matrix, then \code{A = Q T t(Q)}, where \code{Q} is orthogonal, and
\code{T} is upper quasi-triangular (nearly triangular with either 1 by 1 or 2 by 2 
blocks on the diagonal). The eigenvalues of \code{A} are the same as those of \code{T},
which are easy to compute. The Schur form is used most often for computing
non-symmetric eigenvalue decompositions, and for computing functions of
matrices such as matrix exponentials.
}
\references{
Golub, G., and Van Loan, C. F. (1989).
\emph{Matrix Computations,}
2nd edition, Johns Hopkins, Baltimore.
}
\seealso{
\code{\link{schur.Matrix}}
}
\examples{
library(Matrix)
x <- Matrix( rnorm(9), 3, 3)
schur(x) 
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
