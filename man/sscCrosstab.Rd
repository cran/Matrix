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
sscCrosstab(flist, upper = TRUE)
}
\arguments{
  \item{flist}{A list of two or more factors.  All factors must have the
    same length.}
  \item{upper}{logical - if \code{TRUE} the upper triangle of the
    pairwise cross-tabulation matrix is computed and returned.  Otherwise
    the lower triangle is returned.  Defaults to \code{TRUE}.}
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
