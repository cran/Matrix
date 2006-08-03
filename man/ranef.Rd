\name{ranef}
\docType{genericFunction}
\docType{methods}
\alias{ranef}
\alias{ranef-methods}
\alias{ranef,ANY-method}
\alias{ranef,lmer-method}
\alias{ranef,mer-method}
\title{Extract Random Effects}
\usage{
ranef(object, \dots)
\S4method{ranef}{lmer}(object, postVar, \dots)
}
\description{
  A generic function to extract, and optionally accumulate, the random effects.
}
\arguments{
  \item{object}{an object of a class from which random effects
    estimates can be extracted.}
  \item{postVar}{an optional logical argument indicating if the
    conditional variance covariance matrices, also called the
    \dQuote{posterior variances}, of the random effects should be
    included.  Default is \code{FALSE}.}
  \item{\dots}{some methods for this generic function require additional
    arguments.}
}
\value{
  A list of data frames, one for each grouping factor for the random
  effects.  The number of rows in the data frame is the number of levels
  of the grouping factor.  The number of columns is the dimension of the
  random effect associated with each level of the factor.

  If \code{postVar} is \code{TRUE} each of the data frames has an
  attribute called \code{"postVar"} which is a three-dimensional array
  with symmetric faces.
}
\details{
  If grouping factor i has k levels and j random effects per level the ith
  component of the list returned by \code{ranef} is a data frame with k
  rows and j columns.  If \code{postVar} is \code{TRUE} the
  \code{"postVar"} attribute is an array of dimension j by j by k.  The
  kth face of this array is a positive definite symmetric j by j
  matrix.  If there is only one grouping factor in the model the
  variance-covariance matrix for the entire random effects vector,
  conditional on the estimates of the model parameters and on the data
  will be block diagonal and this j by j matrix is the kth diagonal block.
  With multiple grouping factors the faces of the \code{"postVar"}
  attributes are still the diagonal blocks of this conditional
  variance-covariance matrix but the matrix itself is no longer block
  diagonal.
}
\note{
  To produce a \dQuote{caterpillar plot} of the random effects apply
  \code{\link[lattice]{qqmath}} to the result of \code{ranef} with
  \code{postVar = TRUE}.
}
\examples{
data(sleepstudy)
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
fm2 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy)
ranef(fm1)
str(rr1 <- ranef(fm1, postVar = TRUE))
qqmath(rr1)
str(ranef(fm2, postVar = TRUE))
}
\keyword{models}
\keyword{methods}
