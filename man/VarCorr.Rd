\name{VarCorr}
\docType{genericFunction}
\alias{VarCorr}
\alias{VarCorr,mer-method}
\title{Extract variance and correlation components}
\usage{
VarCorr(x, \dots)
}
\arguments{
  \item{x}{a fitted model object, usually an object inheriting from
    class \code{lmer}.
  }
  \item{\dots}{Additional, optional arguments for some methods.  At
    present none are used.}
}
\description{
  Extract the estimated variances, standard
  deviations, and correlations of the random-effects terms in a
  linear mixed-effects model, of class \code{lme}, or a generalized
  linear mixed-effects model. When appropriate,
  the within-group error variance and standard deviation are also
  calculated.
}
\value{
  an object of class \code{VarCorr}.
}
\seealso{\code{\link{lmer}}, \code{\link{VarCorr-class}}}
\examples{
(fm2 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy))
VarCorr(fm2)
}
\keyword{models}
