\name{lmer-class}
\title{Mixed Model Representations and *mer Methods}
\docType{class}
% Classes
\alias{lmer-class}
\alias{mer-class}
\alias{coef.lmer-class}
\alias{summary.mer-class}
\alias{summary.lmer-class}
\alias{ranef.lmer-class}
% Methods
\alias{anova,mer-method}
\alias{coef,mer-method}
\alias{coef,mer-method}
\alias{coerce,lmer,dsTMatrix-method}
\alias{coerce,mer,dtCMatrix-method}
\alias{coerce,mer,pMatrix-method}
\alias{chol,lmer,ANY-method}
\alias{confint,mer-method}
\alias{deviance,mer-method}
\alias{fitted,mer-method}
\alias{formula,mer-method}
\alias{logLik,mer-method}
\alias{plot,coef.lmer-method}
\alias{plot,ranef.lmer-method}
\alias{resid,mer-method}
\alias{residuals,mer-method}
\alias{show,mer-method}
\alias{simulate,mer-method}
\alias{summary,mer-method}
\alias{terms,lmer-method}
\alias{update,mer-method}
\alias{vcov,mer-method}
\alias{with,lmer-method}
%
\alias{deviance,summary.mer-method}
\alias{logLik,summary.mer-method}
\alias{vcov,summary.mer-method}
\alias{summary,summary.mer-method}
\description{
  The \code{mer} class is mixed-effects representation using a
  blocked, sparse, symmetric, column-oriented matrix and dense
  matrices.  The \code{lmer} class extends \code{mer} to represent a
  fitted linear mixed-effects model.  The \code{summary.mer} and
  \code{summary.lmer} classes represent summaries of these objects.  The
  \code{coef.lmer} and \code{ranef.lmer} classes represent fixed and
  random effect \dQuote{coefficients}.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("lmer", ...)}
  or, more commonly via the \code{\link{lmer}} function.%
  %% Doug, will we want to keep this?
  %% or by applying the \code{summary} generic to an \code{lmer} object.
}

\section{Slots}{
  \describe{
    \item{\code{flist}:}{The list of grouping factors for the random
      effects.}
    \item{\code{Zt}:}{The transpose of the model matrix for the random
      effects stored as a \code{\link[dgCMatrix-class]{dgCMatrix}}.}
    \item{\code{X}:}{Model matrix for the fixed effects.}
    \item{\code{y}:}{Numeric response vector.}
    \item{\code{wts}:}{Numeric weights vector.}
    \item{\code{wrkres}:}{Numeric vector of working residuals.}
    \item{\code{method}:}{Character vector of length 1 containing
      \code{"ML"}, for maximum likelihood, or \code{"REML"}, for
      restricted maximum likelihood, or \code{"PQL"}, for penalized
      quasi-likelihood, or \code{"Laplace"}, for the Laplacian
      approximation, or \code{"AGQ"}, for adaptive Gauss-Hermite
      quadrature.}
    \item{\code{useScale}:}{logical indicating if variance-covariance
      matrices include a scale factor.  This is \code{FALSE} for the
      glmm families \code{binomial} and \code{poisson}, otherwise
      \code{TRUE}.}
    \item{\code{family}:}{a copy of the glm family.}
    \item{\code{call}:}{a copy of the call to the model-fitting
	function.}
    \item{\code{cnames}:}{a list of character vectors of column names
      of the random-effects variance-covariance matrix associated with
      each grouping factor and the fixed-effects model matrix.}
    \item{\code{nc}:}{integer vector giving the number of columns the
      variance-covariance matrix associated with each grouping factor.}
    \item{\code{Gp}:}{integer vector of group pointers within the random
      effects and the rows of the transposed model matrix in the
      \code{Zt} slot.  The elements of \code{Gp} are the 0-based index
      of the first element corresponding to each grouping factor.  Thus
      the first element is always 0.  The last element is the total
      length of the random effects vector (also the total number of rows
      in the matrix in the \code{Zt} slot).}
    \item{\code{XtX}:}{The crossproduct matrix of the model matrix in
      the \code{X} slot, stored as a
      \code{\link{dpoMatrix-class}{"dpoMatrix"}} object.}
    \item{\code{ZtZ}:}{The crossproduct matrix of the model matrix for
      the random effects, stored as a \code{\link{dpoMatrix-class}{"dpoMatrix"}} object.
      representation of the original Z'Z matrix.}
    \item{\code{ZtX}:}{The original Z'X matrix, stored as a
      \code{\link{dgeMatrix-class}{"dgeMatrix"}} object.}
    \item{\code{Zty}:}{The numeric vector Z'y.}
    \item{\code{Xty}:}{The numeric vector X'y.}
    \item{\code{Omega}:}{A list of positive-definite matrices stored as
      \code{\link{dpoMatrix-class}{"dpoMatrix"}} objects that are the
      relative precision matrices of the random effects associated with
      each of the grouping factors.}
    \item{\code{L}:}{The sparse, supernodal lower Cholesky factor of
      Z'Z+Omega, stored as a \code{\link[CHMfactor-class]{"dCHMsuper"}} object.}
    \item{\code{RZX}:}{An rectangular matrix, stored as a
      \code{\link[dgeMatrix-class]{"dgeMatrix"}} object, which is the RZX
	component of the Cholesky factor of \code{Z'Z+Omega} at the
	current values of the parameters.}
    \item{\code{RXX}:}{An upper triangular matrix, stored as a
      \code{\link[dtrMatrix-class]{"dtrMatrix"}} object, which is the RXX
	component of the Cholesky factor of \code{Z'Z+Omega} at the
	current values of the parameters.}
    \item{\code{rZy}:}{A numeric vector - the rZy
	component of the Cholesky factor of \code{Z'Z+Omega} at the
	current values of the parameters.}
    \item{\code{rXy}:}{A numeric vector - the rXy
	component of the Cholesky factor of \code{Z'Z+Omega} at the
	current values of the parameters.}
    \item{\code{devComp}:}{Numeric vector of components from which the
      deviance can be calculated.  The named components are: \code{"n"},
      the number of observations (length of the \code{"y"} slot),
      \code{"p"}, the length of the fixed-effects vector (number of
      columns of the matrix in the \code{"X"} slot), \code{"yty"}, the
      sum of squares of the \code{y},\code{"logryy2"}, the logarithm of
      the penalized residual sum of squares, \code{"logDetL2"}, the
      logarithm of the square of the determinant of \code{L},
      \code{"logDetOmega"}, the logarithm of the determinant of
      \code{Omega}, and \code{"logDetRXX"}, the logarithm of the
      determinant of the triangular matrix in the \code{"RXX"} slot.
      This slot's contents are current if \code{status$factored} is \code{TRUE}.}
    \item{\code{deviance}:}{Numeric vector of length 2 containing the
      deviance corresponding to the maximum likelihood (ML) and REML
      criteria.  This slot's contents are current if
      \code{status$factored} is \code{TRUE}.}
    \item{\code{fixef}:}{Numeric vector containing the conditional estimates
      of the fixed effects.  This slot's contents are current if
      \code{status$secondary} is \code{TRUE}.}
    \item{\code{ranef}:}{Numeric vector containing the conditional modes
      of the random effects.  This slot's contents are current if
      \code{status$secondary} is \code{TRUE}.}
    \item{\code{RZXinv}:}{An rectangular matrix, stored as a
      \code{\link[dgeMatrix-class]{"dgeMatrix"}} object, which is the RZX
	component of the inverse of the Cholesky factor of
	\code{Z'Z+Omega} at the current values of the parameters. This
	slot's contents are current if \code{status$gradComp} is
	\code{TRUE}.}
    \item{\code{bVar}:}{A list of the diagonal inner blocks (upper
      triangles only) of the positive-definite matrices on the diagonal
      of the inverse of \code{ZtZ+Omega}. With the appropriate scale
      factor (and conversion to a symmetric matrix) these are the
      conditional variance-covariance matrices of the random effects.}
    \item{\code{gradComp}:}{A list of three-dimensional arrays each with four
      faces.  See the code for the documentation of these components.}
    \item{\code{status}:}{Logical vector of length 4 indicating if the
      object has been factored, if the secondary slots have been
      evaluated, if the gradient components have been evaluated and if
      the Hessian components have been evaluated.}
    \item{\code{call}:}{A copy of the matched call to the function that
      created the object.}
    \item{\code{useScale}:}{logical indicator of the use of the scale
      factor when calculating variance-covariance objects.  Only present
      in the \code{summary.lmer} class.}
    \item{\code{showCorrelation}:}{Logical indicator of whether to
      display the correlation of the fixed-effects estimates when
      showing the object. Only present in the \code{summary.lmer} class.}
  }
  The \code{"summary.mer"} and \code{"summary.lmer"} classes \emph{contain}
  the \code{"mer"} and \code{"lmer"} classes respectively and have
  additional slots,
  \describe{
    \item{\code{isG}:}{logical indicating if we have a
      \bold{g}eneralized linear mixed effects model.}
    \item{\code{methTitle}:}{character string specifying a method title}
    \item{\code{logLik}:}{the same as \code{\link{logLik}(object)}.}
    \item{\code{ngrps}:}{...}% FIXME
    \item{\code{sigma}:}{the scale ...}% FIXME
    \item{\code{coefs}:}{the coefficients matrix...}% FIXME
    \item{\code{vcov}:}{the same as \code{\link{vcov}(object)}.}
    \item{\code{REmat}:}{the formatted Random-Effects matrix}
    \item{\code{AICtab}:}{...}% FIXME
  }% summary slots
}
\usage{
\S4method{logLik}{mer}(object, REML = object@method == "REML", \dots)
\S4method{vcov}{mer}(object, REML = object@method == "REML",
       useScale = object@useScale, \dots)
}
\arguments{
  \item{object}{object of class \code{mer} or \code{lmer}.}
  \item{REML}{logical indicating if REML should be used.}
  \item{useScale}{logical indicating if a scale should be computed and
    used; otherwise \code{sc <- 1} is used.}
  \item{\dots}{potential further arguments passed to methods.}
}
\section{Methods}{
% TODO:
%  "show": object = "lmer" (inherited from object = "mer")
%  "with": data = "lmer"
  \describe{
    \item{anova}{\code{signature(object = "mer")}: model comparisons for
      \code{mer} objects, see \code{\link{anova}}.}
      %% FIXME: say a bit more on lmer - model comparison ??
    \item{coef}{\code{signature(object = "lmer")}: Extract the
      parameters that determine the Omega matrices.  The optional
      logical argument \code{unconstr} determines the constrained or
      unconstrained parameterization.}
    \item{coef<-}{\code{signature(object = "lmer", value = "numeric")}: Assign the
      parameters that determine the Omega matrices.  The optional
      logical argument \code{unconstr} determines the constrained or
      unconstrained parameterization. }
    \item{coerce}{\code{signature(from = "lmer", to = "dtTMatrix")}:
      Ensure that the \code{lmer} object has been factored but not
      inverted and extract the unit lower triangular matrix \eqn{L}
      from the \eqn{LDL} decomposition of
      \eqn{Z^\prime Z+\Omega}{crossprod(Z)+Omega}.}
    \item{coerce}{\code{signature(from = "lmer", to = "dsTMatrix")}:
      Ensure that the \code{lmer} object has been factored but not
      inverted and extract the symmetric matrix
      \eqn{Z^\prime Z+\Omega}{crossprod(Z)+Omega}.}
    \item{deviance}{\code{signature(object = "mer")}: extracts the
      \code{\link{deviance}} of the fitted model.}
%%% some of the following are also in \usage{..} because they have extra args:
    \item{logLik}{\code{signature(object = "mer")}: Extract the
      log-likelihood from of the fitted model, see also \code{\link{logLik}}.}
    \item{update}{\code{signature(object = "mer")}: see
      \code{\link{update}} on how to update fitted models.}
    \item{simulate}{\code{signature(object = "mer")}: see
      \code{\link{simulate}} on how to simulate from fitted models.}
    \item{solve}{\code{signature(a = "lmer", b = "missing")}: Invert
      the decomposed matrices.}
    \item{terms}{\code{signature(x = "lmer")}: Extract
      \code{\link{terms}} component.}
    \item{vcov}{\code{signature(object = "mer")}: Calculate
      variance-covariance matrix of the \emph{fixed} effect terms,
      see also \code{\link{vcov}}.
    }
  }
}
\seealso{
  \code{\link{lmer}()} produces these objects.
  \cr
  \code{\link{VarCorr}} for extracting the variance and
  correlation components of the \emph{random}-effects terms.
}
\examples{
(fm2 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject),
             data = sleepstudy))
logLik(fm2)
(V2 <- vcov(fm2))
}
\keyword{classes}
