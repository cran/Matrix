\name{lmer-class}
\docType{class}
\alias{lmer-class}
\alias{chol,lmer,ANY-method}
\alias{coef,lmer-method}
\alias{coef<-,lmer,numeric-method}
\alias{solve,lmer,missing-method}
\title{lme models representation}
\description{A representation of a linear mixed-effects model using a
  combination of sparse, symmetric, column-oriented matrices and dense matrices.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("lmer", ...)}
or, more commonly, via the \code{lme1} function.}
}
\section{Slots}{
  \describe{
    \item{\code{Omega}:}{A list of numeric matrices providing the
      components of symmetric, positive-definite matrix Omega.  Only the
      upper triangle of each component is used and stored.}
    \item{\code{D}:}{A list of the diagonal factors (lower triangle) in the
      LDL' decomposition of Z'Z+W.}
    \item{\code{T}:}{A list of objects of class \code{cscMatrix}
      providing the pairwise cross-tabulation of the grouping factors.}
    \item{\code{ZZx}:}{A list of arrays comprising ZtZ}
    \item{\code{RXX}:}{A matrix which is the (augmented) RXX component
      or the corresponding component from the inverse of Z'Z+W}
    \item{\code{RZX}:}{A matrix which is the (augmented) RZX component
      or the corresponding component from the inverse of Z'Z+W}
    \item{\code{XtX}:}{The original X'X matrix}
    \item{\code{ZtX}:}{The original Z'X matrix}
    \item{\code{cnames}:}{Column names of the model matrices.}
    \item{\code{deviance}:}{Numeric vector of length 2 containing the
      deviance corresponding to the maximum likelihood (ML) and REML
      criteria.  This slot's contents are current if
      \code{status$factored} is \code{TRUE}.}
    \item{\code{devComp}:}{Numeric vector of length 4 giving the
      components used to calculate the deviance. This slot's contents
      are current if \code{status$factored} is \code{TRUE}.}
    \item{\code{levels}:}{List of the levels of the grouping factors}
    \item{\code{nc}:}{Integer vector with the number of columns in
      (augmented) model matrices and the number of observations in the
      model frame.}
    \item{\code{status}:}{Logical vector of length 2 indicating if the
      object has been factored and if the factorization has been
      inverted.}
  }
}
\section{Methods}{
  \describe{
    \item{chol}{\code{signature(x = "lmer", pivot = "ANY")}: Factor
      the matrix Z'Z+W}
    \item{coef}{\code{signature(object = "lmer")}: Extract the
      parameters that determine the Omega matrices.  The optional
      logical argument \code{unconstr} determines the constrained or
      unconstrained parameterization.}
    \item{coef<-}{\code{signature(object = "lmer", value = "numeric")}: Assign the
      parameters that determine the Omega matrices.  The optional
      logical argument \code{unconstr} determines the constrained or
      unconstrained parameterization. }
    \item{solve}{\code{signature(a = "lmer", b = "missing")}: Invert
      the decomposed matrices.}
  }
}
%\references{}
%\author{}
%\note{}
\seealso{
  \code{\link{sscCrosstab-class}}, \code{\link{sscMatrix-class}},
  \code{\link{tscMatrix-class}}
}
%\examples{}
\keyword{classes}
