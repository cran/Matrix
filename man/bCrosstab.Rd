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
  A list with two components,
  \item{\code{flist}}{a possibly reordered copy of the
    first argument, and}
  \item{\code{ctab}}, itself a list of all pairwise
  crosstabulation matrices.  These include the (diagonal) matrices of
  \dQuote{interactions of a factor with itself}, see the example.
}
%\references{}
%\author{}
%\note{}
\seealso{\code{\link{:}} for factor interaction,
  \code{\link{xtabs}} for another cross tabulation function.
}
\examples{
str(esoph)
bCrosstab(with(esoph,list(agegp, alcgp))) $ ctab

grps <- with(esoph, list(age = agegp, alc = alcgp, tob = tobgp))
bCn <- bCrosstab( grps, reorder = FALSE)
stopifnot(all.equal(bCn$flist, grps))
str(bC <- bCrosstab(grps))
bC$ctab$"alc:age" # generally 4 obs. per cell, but not always:
aa.m <- as(bC$ctab$"alc:age", "matrix")
table(aa.m) # #(2,3,4) = (2, 4, 18)
stopifnot(all.equal(as.vector(t(aa.m)),
                    as.vector(table(with(esoph, alcgp:agegp))), tol=0))
}
\keyword{manip}
\keyword{algebra}
