\name{schurmod}
\alias{schurmod}
\title{
Reordered Schur Factorization
}
\description{
Reorders the Schur factorization of a square matrix.
}
\usage{
schurmod(x, indexes=<see below>, vectors=<see below>)
}
\arguments{
\item{x}{
an object of class \code{schur.Matrix}.
}
\item{indexes}{
a integer vector of length 2 indicating which indexes of the eigenvalues 
involved in the reordering. If \code{indexes = c(i,j)}, then the block or element
in row \code{i} is moved to row \code{j} by a sequence of transpositions between
adjacents blocks or elements.
}
\item{vectors}{
a logical value telling whether or not to update and return the Schur vectors
if they are present in \code{x}. The default is to update and return the Schur 
vectors.
}}
\value{
An object of class \code{schur.Matrix} containing the reordered Schur form and 
(optionally) the corresponding Schur vectors.
}
\details{
Based on the functions dtrexc and ztrexc from Lapack (Anderson et al. 1994).
}
\references{
Anderson, E., et al. (1994).
\emph{LAPACK User's Guide,}
2nd edition, SIAM, Philadelphia.
}
\seealso{
\code{\link{schur.Matrix}}, \code{\link{schur.Matrix.object}}
}
\section{BUGS}{
Does not return the proper ordering in all cases. A bug report will be filed
with the Lapack developers.
}
\examples{
n <- 9
x <- Matrix( rnorm( n*n, sd = 100), nrow = n, ncol = n) 
schur.x <- schur(x)                                 # schur factorization
schur.x <- schurmod(schur.x, indexes = c(1,n))      # reordered factorization
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
