\name{bCrosstab}
\alias{bCrosstab}
\title{Create pairwise crosstabulation}
\description{
  Create the pairwise crosstabulation of two or more factors.
}
\usage{
bCrosstab(flist, reorder = TRUE)
}
\arguments{
  \item{flist}{A list of two or more factors.  All factors must have the
    same length.}
  \item{reorder}{Logical flag determining if reordering of the factors
    is allowed.}
}
%\details{}
\value{
  A list with components \code{flist}, a possibly reordered copy of the
  first argument, and \code{ctab}, the pairwise crosstabulation.
}
%\references{}
%\author{}
%\note{}
%\seealso{}
%\examples{}
\keyword{manip}
\keyword{algebra}
