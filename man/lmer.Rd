\name{lmer}
\docType{genericFunction}
\alias{lmer}
\alias{lmer,formula-method}
\concept{GLMM}% << so it's found
\title{Fit (Generalized) Linear Mixed-Effects Models}
\description{
  This generic function fits a linear mixed-effects model with nested or
  crossed grouping factors for the random effects.
}
\synopsis{
lmer(formula, data, family, method = c("REML", "ML", "PQL", "Laplace", "AGQ"),
    control = list(), subset, weights, na.action, offset,
    model = TRUE, x = FALSE, y = FALSE, \dots)
}
\usage{
lmer(formula, data, family,
     method = c("REML", "ML", "PQL", "Laplace", "AGQ"),\dots)
}
\arguments{
  \item{formula}{a two-sided linear formula object describing the
    fixed-effects part of the model, with the response on the left of a
    \code{~} operator and the terms, separated by \code{+} operators, on
    the right.  The vertical bar character \code{"|"} separates an
    expression for a model matrix and a grouping factor.}
  \item{data}{an optional data frame containing the variables named in
    \code{formula}.  By default the variables are taken from the
    environment from which \code{lmer} is called.}
  \item{family}{
    a GLM family, see \code{\link[stats]{glm}}.  If \code{family} is
    missing then a linear mixed model is fit; otherwise a generalized
    linear mixed model is fit.}
  \item{method}{a character string.  For a linear mixed model the
    default is \code{"REML"} indicating that the model should be fit by
    maximizing the restricted log-likelihood.  The alternative is
    \code{"ML"} indicating that the log-likelihood should be
    maximized.  (This method is sometimes called "full" maximum
    likelihood.)  For a generalized linear mixed model the criterion
    is always the log-likelihood but this criterion does not have a
    closed form expression and must be approximated.  The default approximation is
    \code{"PQL"} or penalized quasi-likelihood.  Alternatives are
    \code{"Laplace"} or \code{"AGQ"} indicating the Laplacian and
    adaptive Gaussian quadrature approximations respectively.  The
    \code{"PQL"} method is fastest but least accurate.  The
    \code{"Laplace"} method is intermediate in speed and accuracy.
    The \code{"AGQ"} method is the most accurate but can be
    considerably slower than the others.}
  \item{\dots}{Optional arguments for methods.  Currently none are used.}
  }
}
\details{
  This is a revised version of the \code{\link[nlme]{lme}} function from the
  \pkg{nlme} package.  This version uses a different method of
  specifying random-effects terms and allows for fitting generalized
  linear mixed models as well as linear mixed models.

  Additional standard arguments to model-fitting functions can be passed
  to \code{lmer}.
  \describe{
    \item{subset}{an optional expression indicating the subset of the rows of
      \code{data} that should be used in the fit. This can be a logical
      vector, or a numeric vector indicating which observation numbers are
      to be included, or a  character  vector of the row names to be
      included.  All observations are included by default.}
    \item{na.action}{a function that indicates what should happen when the
      data contain \code{NA}s.  The default action (\code{na.fail}) causes
      \code{lme} to print an error message and terminate if there are any
      incomplete observations.}
    \item{control}{a named list of control values for the estimation
      algorithm.  Possible control options and their values are:
      \describe{
	\item{maxIter}{maximum number of iterations for the \code{lme}
	  optimization algorithm. Default is 50.}
	\item{msMaxIter}{maximum number of iterations
	  for the \code{nlminb} optimization step inside the \code{lme}
	  optimization. Default is 200.}
	\item{niterEM}{number of iterations for the EM algorithm used to refine
	  the initial estimates of the random effects variance-covariance
	  coefficients. Default is 15.}
	\item{msVerbose}{a logical value passed as the \code{trace} argument to
	  \code{nlminb} (see documentation on that function). Default is
	  \code{getOption("verbose")}.}
	\item{PQLmaxIt}{maximum number of iterations for the PQL
	  algorithm when fitting generalized linear mixed models.
	  Default is 30.}
	\item{EMverbose}{a logical value indicating if verbose output
	  should be produced during the EM iterations.  Default is
	  \code{getOption("verbose"}.}
	\item{analyticGradient}{a logical value indicating if the
	  analytic gradient of the objective should be used.  This
	  option is for testing purposes and would not normally be
	  changed from the default.  Default is \code{TRUE}.}
	\item{analyticHessian}{a logical value indicating if the
	  analytic Hessian of the objective should be calculated. This
	  is an experimental feature and at present the default is
	  \code{FALSE}. In the future we may use the analytic Hessian in
	  the optimization.}
      }
      Defaults to an empty list.}
    \item{model, x}{logicals.  If \code{TRUE} the corresponding
      components of the fit (the model frame, the model matrices)
      are returned.}
  }
}
\value{
  An \code{\link[Matrix]{lmer-class}{lmer}} object.
}
%\references{ ~put references to the literature/web site here ~ }
%\author{}
%\note{}
\seealso{\code{\link[Matrix]{lmer-class}}, \code{\link[stats]{lm}}}
\examples{
(fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy))
(fm2 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy))
anova(fm1, fm2)
}
\keyword{models}
\keyword{methods}
