\name{ssclme-class}
\docType{class}
\alias{ssclme-class}
\alias{coerce,ssclme,sscMatrix-method}
\alias{coerce,ssclme,tscMatrix-method}
\title{lme models as sparse, symmetric matrices}
\description{A representation of a linear mixed-effects model using a
  combination of sparse, symmetric, column-oriented matrices and dense matrices.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ssclme", ...)}
or, more commonly, via the \code{ssclme} constructor function.}
}
\section{Slots}{
  \describe{
    \item{\code{Gp}:}{integer vector of pointers to the beginning of
      each group of columns.  Its length is 1 greater than the number of
      groups.  Successive differences give the number of levels of each
      factor.} 
    \item{\code{p}:}{as for \code{"sscMatrix"} objects.}
    \item{\code{i}:}{as for \code{"sscMatrix"} objects.}
    \item{\code{x}:}{as for \code{"sscMatrix"} objects.}
    \item{\code{Dim}:}{as for \code{"sscMatrix"} objects.}
    \item{\code{nc}:}{An integer vector of the number of columns in each
      model matrix.  The last entry is \code{p+1} where \code{p} is the
      number of columns in the fixed-effects model matrix.}
    \item{\code{ZtX}:}{A (dense) matrix equivalent to \code{crossprod(Z,X)}.}
    \item{\code{XtX}:}{A (dense) matrix equivalent to the upper triangle
      of \code{crossprod(X)} where \code{X} has the response column appended.}
    \item{\code{Lp}:}{An integer vector of pointers to the beginnings of
      the columns of the unit lower triangular factor L in the
      LDL' decomposition of Z'Z+W.}
    \item{\code{Parent}:}{An integer vector from the symbolic analysis
      of the cross-tabulation.  It is used to compute L efficiently.}
    \item{\code{Li}:}{Row indices of L.}
    \item{\code{D}:}{Diagonal of the diagonal matrix in the
      LDL' decomposition of Z'Z+W.}
    \item{\code{DIsqrt}:}{Square root of the inverse of the diagonal of
      D in LDL'.}
    \item{\code{Lx}:}{Non-zero elements in the (strict) lower triangle
      of L.}
    \item{\code{RZX}:}{A dense matrix corresponding to \code{ZtX} in the
      decomposition.}
    \item{\code{RXX}:}{A dense matrix corresponding to \code{XtX} in the
      decomposition.}
    \item{\code{LIp}:}{An integer vector of pointers to
      the columns of the unit lower triangular inverse of L.}
    \item{\code{LIi}:}{Row indices of L-inverse.}
    \item{\code{LIx}:}{Non-zero elements in the (strict) lower triangle
      of L-inverse.}
    \item{\code{status}:}{Logical vector of length 2 indicating if the
      object has been factored and if the factorization has been
      inverted.}
    \item{\code{deviance}:}{Numeric vector of length 2 containing the
      deviance corresponding to the maximum likelihood (ML) and REML
      criteria.  This slot's contents are current if
      \code{status$factored} is \code{TRUE}.}
    \item{\code{devComp}:}{Numeric vector of length 4 giving the
      components used to calculate the deviance. This slot's contents
      are current if \code{status$factored} is \code{TRUE}.}
    \item{\code{bVar}:}{A list of numeric 3-dimensional arrays giving
      the conditional variances of the random effects. This slot's
      contents are current if \code{status$inverted} is \code{TRUE}.}
    \item{\code{Omega}:}{A list of numeric matrices providing the
      components of symmetric, positive-definite matrix Omega.  Only the
      upper triangle of each component is used and stored.}
  }
}
\section{Methods}{
No methods defined with class "ssclme" in the signature.
}
%\section{Methods}{}
%\references{}
%\author{}
%\note{}
\seealso{
  \code{\link{sscCrosstab-class}}%, \code{\link{ssclme}}
}
%\examples{}
\keyword{classes}
