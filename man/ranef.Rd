% $Id: random.effects.Rd,v 1.1 2002/11/27 22:49:45 bates Exp $
\name{ranef}
\docType{genericFunction}
\docType{methods}
\alias{ranef}
\alias{random.effects}
\alias{ranef-methods}
\alias{ranef,ANY-method}
\alias{ranef,lmer-method}
\alias{ranef,mer-method}
\title{Extract Random Effects}
\usage{
ranef(object, \dots)
\S4method{ranef}{lmer}(object, accumulate, \dots)
random.effects(object, \dots)   # deprecated
}
\description{
  A generic function to extract, and optionally accumulate, the random effects.
}
\arguments{
 \item{object}{an object of a class from which random effects
   estimates can be extracted.}
 \item{accumulate}{an optional logical argument indicating if the random
   effects for nested grouping factors should be accumulated.}
 \item{\dots}{some methods for this generic function require additional
   arguments.}
}
\value{
  a list of matrices.
}
\examples{
ranef(lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy))
}
\keyword{models}
\keyword{methods}
