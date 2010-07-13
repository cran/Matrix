\name{respModule-class}
\Rdversion{1.1}
\docType{class}
\alias{respModule-class}
\alias{glmRespMod-class}
\alias{nlsRespMod-class}
\alias{nglmRespMod-class}
\title{"respModule" and derived classes}
\description{
  The \code{"respModule"} class is the virtual base class of
  response modules for \code{\linkS4class{glpModel}} model objects.
  Classes that inherit from \code{"respModule"} include
  \code{\linkS4class{glmRespMod}}, for
  generalized linear models, \code{\linkS4class{nlsRespMod}}, for
  nonlinear models and \code{\linkS4class{nglmRespMod}} for generalized
  nonlinear models.
}
\section{Objects from the Class}{
  Objects from these classes are usually created with
  \code{\link{mkRespMod}} as part of an \code{\linkS4class{glpModel}}
  object returned by model-fitting functions such as the hidden function
  \code{glm4}.
}
\section{Slots}{
  \describe{
    \item{\code{mu}:}{Fitted mean response.}
    \item{\code{offset}:}{offset in the linear predictor -- always
      present even if it is a vector of zeros.  In an
      \code{\linkS4class{nlsRespMod}} object the length of the offset can be
      a multiple of the length of the response.}
    \item{\code{sqrtXwt}:}{the matrix of
    weights for the model matrices, derived from the \code{sqrtrwt} slot.}
    \item{\code{sqrtrwt}:}{Numeric vector of the square roots
      of the weights for the residuals.  For
      \code{respModule} and \code{\linkS4class{nlsRespMod}}
      objects these are constant.  For
      \code{\linkS4class{glmRespMod}} and \code{\linkS4class{nglmRespMod}}
      objects these are updated at each iteration of the iteratively
      reweighted least squares algorithm.}
    \item{\code{weights}:}{Prior
      weights -- always present even when it is a vector of ones.}
    \item{\code{y}:}{Numeric response vector.}
    \item{\code{family}:}{a glm family, see \code{\link{family}} for
      details - \code{glmRespMod} objects only.}
    \item{\code{eta}:}{numeric vector, the linear predictor that is
      transformed to the conditional mean via the link function -
      \code{glmRespMod} objects only.}
    \item{\code{n}:}{a numeric vector used for calculation of the aic family
      function (it is really only used with the binomial family but
      we need to include it everywhere) - \code{glmRespMod} objects only.}
    \item{\code{nlenv}:}{an environment in which to evaluate the
      nonlinear model function - \code{nlsRespMod} objects only.}
    \item{\code{nlmod}:}{an unevaluated call to the
      nonlinear model function - \code{nlsRespMod} objects only.}
    \item{\code{pnames}:}{a character vector of parameter names -
      \code{nlsRespMod} objects only.}
  }
}
\section{Methods}{
No methods defined with class "respModule" in the signature.
}
\seealso{\code{\link{mkRespMod}}}
\examples{
showClass("respModule")
showClass("glmRespMod")
showClass("nlsRespMod")
}
\keyword{classes}
