\name{hatTrace}
\alias{hatTrace}
\title{Extract the trace of the `hat' matrix}
\usage{
hatTrace(x)
}
\description{
  Extract the trace of the `hat' matrix from a fitted linear mixed model
  or generalized linear mixed model.
}
\arguments{
 \item{x}{a fitted model object that inherits from the
   \code{\linkS4class{mer}} class}
}
\value{
  The trace of the hat matrix at the parameter estimates.  This is a
  measure of the complexity of the model, conditional on the relative
  variance components.
}
\references{Spiegelhalter, D.J., Best, N.G., Carlin, B.P. and van der
  Linde, A. (2002), "Bayesian measures of model complexity and fit (with
  discussion)", \emph{Journal of the Royal Statistical Society, part B},
  \bold{64}, (4), 583-639.}
\examples{
data(sleepstudy)
hatTrace(lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy))
}
\keyword{models}
\keyword{methods}
