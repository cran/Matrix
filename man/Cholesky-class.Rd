\name{Cholesky-class}
\docType{class}
\alias{Cholesky-class}

\title{Cholesky decompositions}
\description{The \code{"Cholesky"} class is the class of Cholesky
  decomposition of a positive-semidefinite, real matrices}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("Cholesky", ...)}
  or by calls of the form \code{chol(pm)} where \code{pm} inherits from
  the \code{"poMatrix"} class or as a side-effect of other functions
  applied to \code{"poMatrix"} objects.
}
\section{Slots}{
  \describe{
    \item{\code{uplo}:}{inherited from the \code{"trMatrix"} class.}
    \item{\code{diag}:}{inherited from the \code{"trMatrix"} class.}
    \item{\code{x}:}{inherited from the \code{"trMatrix"} class.}
    \item{\code{Dim}:}{inherited from the \code{"trMatrix"} class.}
    \item{\code{rcond}:}{inherited from the \code{"trMatrix"} class.}
    \item{\code{factorization}:}{inherited from the \code{"trMatrix"}
      class.  For this class the \code{factorization} slot is always an
      empty list.}
  }
}
\section{Extends}{
Class \code{"trMatrix"}, directly.
Class \code{"geMatrix"}, by class \code{"trMatrix"}.
Class \code{"Matrix"}, by class \code{"trMatrix"}.
}
\section{Methods}{
No methods defined with class "Cholesky" in the signature.
}
\seealso{
  \code{\link{trMatrix-class}}, \code{chol}, \code{poMatrix-class}
}
\examples{
print(pm <- crossprod(Matrix(rnorm(18), nrow = 6, ncol = 3)))
print(ch <- chol(pm))
if (toupper(ch@uplo) == "U") crossprod(ch)
}
\keyword{classes}
\keyword{algebra}
