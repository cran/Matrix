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
\usage{
lmer(formula, data, family, method,
     control, start, subset, weights, na.action, offset,
     model, x, y, \dots)

\S4method{lmer}{formula}(formula, data, family, method,
    control, start, subset, weights, na.action, offset,
    model, x, y, \dots)
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
  \item{control}{a list of control parameters. See below for details.}
  \item{start}{a list of relative precision matrices for the random
    effects.  This has the same form as the slot \code{"Omega"} in a
    fitted model.  Only the upper triangle of these symmetric matrices
    should be stored.}
  \item{subset, weights, na.action, offset}{further model specification
    arguments as in \code{\link[stats]{lm}}; see there for details.}
  \item{model, x, y}{logicals indicating if the corresponding model component
    should be returned.  Note that \code{x} and \code{y} are not yet
    implemented.}
  \item{\dots}{potentially further arguments for methods.  Currently
    none are used.}
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
    \item{control}{a named list of control parameters for the estimation
      algorithm, specifying only the ones to be changed from their
      default values.  Hence defaults to an empty list.\cr
      Possible control options and their default values are:
      \describe{
	\item{\code{maxIter}:}{maximum number of iterations for the \code{lme}
	  optimization algorithm.  Default is 50.}
	\item{\code{tolerance}:}{relative convergence tolerance for the
	  \code{lme} optimization algorithm.  Default is
	  \code{sqrt(.Machine$double.eps)}.}
	\item{\code{msMaxIter}:}{maximum number of iterations
	  for the \code{nlminb} optimization step inside the \code{lme}
	  optimization.  Default is 200.}
	\item{\code{msVerbose}:}{a logical value passed as the
	  \code{trace} argument to \code{nlminb} (see documentation on
	  that function).  Default is \code{getOption("verbose")}.}
	\item{\code{niterEM}:}{number of iterations for the EM algorithm
	  used to refine the initial estimates of the random effects
	  variance-covariance coefficients.  Default is 15.}
	\item{\code{EMverbose}:}{a logical value indicating if verbose output
	  should be produced during the EM iterations.  Default is
	  \code{getOption("verbose")}.}
	\item{\code{PQLmaxIt}:}{maximum number of iterations for the PQL
	  algorithm when fitting generalized linear mixed models.
	  Default is 30.}%% unused currently?
	\item{\code{gradient}:}{a logical value indicating if the
	  analytic gradient of the objective should be used.  Use of an
	  analytic gradient results in more stable convergence but can take
	  longer.  For models with multiple grouping factors the
	  difference in time can be substantial.
	  Default is \code{TRUE}.}
	\item{\code{Hessian}:}{a logical value indicating if the
	  analytic Hessian of the objective should be calculated.  This
	  is an experimental feature and at present the default is
	  \code{FALSE}.  In the future we may use the analytic Hessian in
	  the optimization.}
      }
    }
    \item{model, x}{logicals.  If \code{TRUE} the corresponding
      components of the fit (the model frame, the model matrices)
      are returned.}
  }
}
\value{
  An object of class \code{"\link[Matrix:lmer-class]{lmer}"}.
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
