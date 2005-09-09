\name{simulate}
\docType{genericFunction}
\alias{simulate}
\alias{simulate,lmer-method}
\title{Simulate responses from a fitted model}
\description{
  This generic function generates a sample of responses from a fitted
  model object using the parameter estimates from the fit.  Typically
  these simulated responses are used in a parametric bootstrap sample of
  a test statistic.
}
\usage{
simulate(object, nsim, seed, \dots)
}
\arguments{
  \item{object}{an object representing a fitted model.}
  \item{nsim}{number of response vectors to simulate.  Defaults to 1.}
  \item{seed}{an integer that will be used in a call to \code{set.seed}
    before simulating the response vectors.  The value is saved as an
    attribute named \code{"seed"} in the returned value.  Defaults to a
    random integer between \code{0} and \code{.Machine$integer.max}.}
  \item{\dots}{additional optional arguments.}
}
\value{
  A data frame with \code{nsim} columns and \code{nobs}, the number
  of observations in the fitted model \code{object}, rows.
}
\section{Methods}{
  \describe{
    \item{object = "lmer"}{simulate responses from the fitted model \code{"object"}.}
  }
}
\examples{
(fm2 <- lmer(Reaction ~ Days + (1|Subject) + (0 + Days|Subject), sleepstudy))
sim2 <- simulate(fm2, nsim = 4)
head(sim2)
tail(sim2)
}
\keyword{methods}
\keyword{datagen}
