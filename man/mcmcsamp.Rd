\name{mcmcsamp}
\docType{genericFunction}
\alias{mcmcsamp}
\alias{mcmcsamp,mer-method}
\title{Generate an MCMC sample}
\description{
  This generic function generates a sample from the posterior
  distribution of the parameters of a fitted model using Markov Chain
  Monte Carlo methods.
}
\usage{
mcmcsamp(object, n, verbose, \dots)
}
\arguments{
  \item{object}{An object of a suitable class - usually an
    \code{\linkS4class{lmer}} object.
  }
  \item{n}{integer - number of samples to generate. Defaults to 1.}
  \item{verbose}{logical - if \code{TRUE} verbose output is printed.
    Defaults to \code{FALSE}.}
  \item{\dots}{Some methods for this generic function may take
    additional, optional arguments.  The method for
    \code{\linkS4class{lmer}} objects takes the optional argument
    \code{saveb} which, if \code{TRUE}, causes the values of the random
    effects in each sample to be saved.  Note that this can result in
    very large objects being saved.  Use with caution. A second optional
    argument is \code{trans} which, if \code{TRUE} (the default), returns
    a sample of transformed parameters.  All variances are expressed on
    the logarithm scale and any covariances are converted to Fisher's "z"
    transformation of the corresponding correlation.}
}
\value{
  An object of (S3) class \code{"mcmc"} suitable for use with the
  functions in the "coda" package.
}
\section{Methods}{
  \describe{
    \item{object = "lmer"}{generate MCMC samples from the posterior
      distribution of the parameters of a linear mixed model or a
      generalized linear mixed model.  The prior on the fixed effects
      parameters is taken to be locally uniform.  The prior on the
      variance-covariance matrices of the random effects is taken to be
      the locally non-informative prior described in Box and Tiao
      (1973). Conditional on the current values of the random effects
      these are sampled from a Wishart distribution.}
  }
}
\examples{
require("lattice", quietly = TRUE, character = TRUE)
(fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy))
set.seed(101)
samp1 <- mcmcsamp(fm1, n = 1000)
frm <-
    data.frame(vals = c(samp1), iter = rep(1:nrow(samp1), ncol(samp1)),
    par = factor(rep(1:ncol(samp1), each = nrow(samp1)),labels = colnames(samp1)))
densityplot(~ vals | par, frm, plot = FALSE,
            scales = list(relation = 'free', x = list(axs='i')))
xyplot(vals ~ iter | par, frm, layout = c(1, ncol(samp1)),
       scales = list(x = list(axs = "i"), y = list(relation = "free")),
       main = "Trace plot", xlab = "Iteration number", ylab = "",
       type = "l")
qqmath(~ vals | par, frm, type = 'l',
       scales = list(y = list(relation = 'free')))
if (require("coda", quietly = TRUE, character = TRUE)) {
   print(summary(samp1))
   print(autocorr.diag(samp1))
}
(eDF <- mean(samp1[,"deviance"]) - deviance(fm1)) # potentially useful approximate D.F.
\dontshow{
stopifnot(abs(eDF - 15) < 1)
}
}
\keyword{methods}
\keyword{datagen}
