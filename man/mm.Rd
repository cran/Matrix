\name{mm}
\alias{mm}
\docType{data}
\title{A sample sparse model matrix}
\description{
  A model matrix used in an example by Koenker and Ng.  This is a sparse
  matrix with 1850 rows and 712 columns but only 8758 non-zero entries.
  It is a \code{"cscMatrix"} object.
}
\usage{data(mm)}
%\details{}
%\source{}
\references{
  Roger Koenker and Pin Ng (2003), ``SparseM: A sparse matrix package
  for R'', \emph{J. of Statistical Software}, 8(6).
}
\examples{
data(mm)
data(y)
class(mm)
dim(mm)
str(mm)
str(y)
sparse.sol = solve(crossprod(mm), crossprod(mm, y))
str(sparse.sol)
}
\keyword{datasets}
