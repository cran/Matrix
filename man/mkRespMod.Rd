\name{mkRespMod}
\alias{mkRespMod}
\title{Create a respModule object}
\description{
  Create a \code{\linkS4class{respModule}} object, which could be from a
  derived class such as \code{\linkS4class{glmRespMod}} or
  \code{\linkS4class{nlsRespMod}}.
}
\usage{
mkRespMod(fr, family = NULL, nlenv = NULL, nlmod = NULL)
}
\arguments{
  \item{fr}{
    a model frame, usually created by a call to \code{\link{model.frame}}.}
  \item{family}{
    an optional glm \code{\link{family}} object
    (\code{\linkS4class{glmRespMod}} objects only).}
  \item{nlenv}{
    an environment for evaluation of the nonlinear model, \code{nlmod}.
    (\code{\linkS4class{nlsRespMod}} objects only).}
  \item{nlmod}{
    the nonlinear model function, as a function call
    (\code{\linkS4class{nlsRespMod}} objects only).}
}
\details{
  The internal representation of a statistical model based on
  a linear predictor expression is derived from a \code{\link{formula}}
  expression and a \code{data} argument, possibly supplemented with a
  \code{\link{family}} object and/or a nonlinear model expression.  The
  steps to obtain this representation usually involve calls to
  \code{\link{model.frame}} and to \code{\link{model.matrix}} or
  \code{\link{model.Matrix}}, which encapsulate important parts of this
  process.   This function encapsulates other operations related to
  weights and offsets and to the model family to create a
  \code{\linkS4class{respModule}} object.
}
\value{
  an object of a class inheriting from \code{\linkS4class{respModule}}.
}
%\references{}
%\author{Doug Bates}
\seealso{
  The \code{\linkS4class{respModule}} class description.
}
\examples{
  ## see  help("glpModel-class")
}
\keyword{models}
