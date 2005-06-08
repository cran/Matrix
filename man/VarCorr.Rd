\name{VarCorr}
\docType{genericFunction}
\alias{VarCorr}
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
(fm1 <- lmer(decrease ~ treatment + (1|rowpos) + (1|colpos),
             OrchardSprays))
VarCorr(fm1)
}
\keyword{models}
