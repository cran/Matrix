\name{pedigree-class}
\docType{class}
\alias{pedigree-class}
\alias{chol,pedigree-method}
\alias{coerce,pedigree,dtCMatrix-method}
\alias{coerce,pedigree,data.frame-method}
\alias{head,pedigree-method}
\alias{show,pedigree-method}
\alias{tail,pedigree-method}
%
\title{Class "pedigree"}
\description{
  Objects of class \code{"pedigree"} represent a set of
  individuals that can have two parents including their parent-child
  relations.  The terminology has been taken from cattle breeding.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("pedigree",
    ...)}, or more conveniently, \code{\link{pedigree}(sire= ., dam = ., label = .)}.
}
\section{Slots}{
  \describe{
    \item{\code{sire}:}{integer vector (with some \code{NA} entries),
      denoting a \emph{previous} entry in the pedigree corresponding to
      the current entry's \dQuote{father}.}
    \item{\code{dam}:}{similarly as \code{sire} for the \dQuote{mother}
      of each entry.}
    \item{\code{label}:}{a \code{"character"} vector of the same length
      as \code{sire} and \code{dam} giving a unique ID for the
      corresponding entry.}
  }
}
\section{Methods}{
  \describe{
    \item{chol}{\code{signature(x = "pedigree"}}: returns the
      lower-triangular, sparse matrix (a
      \code{\linkS4class{dtCMatrix}} object) that is
      the left Cholesky factor of the relationship matrix of the pedigree.}
    \item{coerce}{\code{signature(from = "pedigree", to = "dtCMatrix")}: ... }
    \item{coerce}{\code{signature(from = "pedigree", to = "data.frame")}: ... }
    \item{head}{\code{signature(x = "pedigree")}: ... }
    \item{show}{\code{signature(object = "pedigree")}: ... }
    \item{tail}{\code{signature(x = "pedigree")}: ... }
  }
}
\references{
  R. A. Mrode, \emph{Linear Models for the Prediction of Animal Breeding
    Values, 2nd ed}, CABI Publishing, 2005.
}
\seealso{
  \code{\link{pedigree}}
}
\examples{
## Rather use, pedigree()!   The following is "raw code":
## Example from chapter 2 of Mrode (2005)
p1 <- new("pedigree",
          sire = as.integer(c(NA,NA,1, 1,4,5)),
          dam  = as.integer(c(NA,NA,2,NA,3,2)),
          label = as.character(1:6))
p1
(dtc <- as(p1, "dtCMatrix"))   # T-inverse in Mrode's notation
solve(dtc)
chol(p1)  # (left) Cholesky factor of A, the relationship matrix
}
\keyword{classes}
