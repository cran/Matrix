\name{sscCrosstab}
\alias{sscCrosstab}
\title{Create pairwise crosstabulation}
\description{
  Create the pairwise crosstabulation of two or more factors and store
  it as an object of class \code{"sscCrosstab"} which consists of a
  symmetric, sparse, positive semidefinite matrix, a vector of the
  counts of the levels in each factor, and, optionally, a permutation
  and inverse permutation that will minimize fill-in in the Cholesky
  decomposition of this matrix while maintaining separation of the
  groups of columns corresponding to different factors.
}
\usage{
sscCrosstab(flist, permutations = TRUE)
}
\arguments{
  \item{flist}{A list of two or more factors.  All factors must have the
    same length.}
  \item{permutations}{logical - if \code{TRUE}, the default, then a
    permutation to minimize fill-in in the Cholesky decomposition while
    maintaining separation of the groups, and its inverse, are calculated.}
}
%\details{}
\value{
  An object of class \code{"sscCrosstab"}.
}
%\references{}
%\author{}
%\note{}
\seealso{\code{\link{sscCrosstab-class}}}
%\examples{}
\keyword{manip}
\keyword{algebra}
