\name{ssclme-class}
\docType{class}
\alias{ssclme-class}
\alias{coef,ssclme-method}
\alias{coef<-,ssclme-method}
\alias{coerce,ssclme,dsCMatrix-method}
\alias{coerce,ssclme,dtCMatrix-method}
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
    \item{\code{p}:}{as for \code{"dsCMatrix"} objects.}
    \item{\code{i}:}{as for \code{"dsCMatrix"} objects.}
    \item{\code{x}:}{as for \code{"dsCMatrix"} objects.}
    \item{\code{Dim}:}{as for \code{"dsCMatrix"} objects.}
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
  \describe{
    \item{coef}{\code{signature(object = "ssclme")}: Extract the
      "coefficients" from the object as a numeric vector.  If the
      optional argument \code{unconst} is \code{FALSE} (the default),
      the values in the upper triangles of the \code{Omega} matrices are
      concatenated.  Within each matrix the diagonal is given first then
      the elements above the diagonal are given row-wise.  If
      \code{unconst} is \code{TRUE}, an unconstrained parameterization
      based on the LDL decomposition (logarithms of diagonal elements of
      D followed by elements of the strict lower triangle of L
      column-wise) are given.}
    \item{coef<-}{\code{signature(object = "ssclme", value =
	"numeric")}: Assign the \code{Omega} matrices from a numeric
      vector.  The assignment depends on the optional argument
      \code{unconst} as described above.}
    \item{coerce}{\code{signature(from = "ssclme", to = "dsCMatrix")}:
      return the current Z'Z matrix as an \code{dsCMatrix} object.}
    \item{coerce}{\code{signature(from = "ssclme", to = "dtCMatrix")}:
      return the current L matrix from the LDL decomposition of Z'Z as
      a \code{dtCMatrix} object.}
  }
}
%\references{}
%\author{}
%\note{}
\seealso{
  \code{\link{sscCrosstab-class}}, \code{\link{dsCMatrix-class}},
  \code{\link{dtCMatrix-class}}
}
%\examples{}
\keyword{classes}
