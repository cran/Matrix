\name{Cholesky-class}
\docType{class}
\alias{Cholesky-class}
\title{Cholesky Decompositions}
\description{The \code{"Cholesky"} class is the class of Cholesky
  decompositions of positive-semidefinite, real matrices.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("Cholesky", ...)}
  or by calls of the form \code{chol(pm)} where \code{pm} inherits from
  the \code{"dpoMatrix"} class or as a side-effect of other functions
  applied to \code{"dpoMatrix"} objects (see \code{\link{dpoMatrix-class}}).
}
\section{Slots}{
  A Cholesky decomposition is basically a triangular matrix extending
  the \code{"dtrMatrix"} class.
  \describe{
    \item{\code{uplo}:}{inherited from the \code{"dtrMatrix"} class.}
    \item{\code{diag}:}{inherited from the \code{"dtrMatrix"} class.}
    \item{\code{x}:}{inherited from the \code{"dtrMatrix"} class.}
    \item{\code{Dim}:}{inherited from the \code{"dtrMatrix"} class.}
    \item{\code{Dimnames}:}{inherited from the \code{"dtrMatrix"} class.}
    \item{\code{rcond}:}{inherited from the \code{"dtrMatrix"} class.}
    \item{\code{factors}:}{inherited from the \code{"dtrMatrix"}
      class.  For this class the \code{factors} slot is always an
      empty list.}
  }
}
\section{Extends}{
  Class \code{"dtrMatrix"}, directly.
  Class \code{"dgeMatrix"}, by class \code{"dtrMatrix"}.
  Class \code{"Matrix"}, by class \code{"dtrMatrix"}.
}
\section{Methods}{
  No methods defined with class "Cholesky" in the signature.
}
\seealso{
  \code{\link{dtrMatrix-class}}, \code{\link{chol}},
  \code{\link{dpoMatrix-class}}.
}
\examples{
(pm <- crossprod(Matrix(rnorm(18), nrow = 6, ncol = 3)))
(ch <- chol(pm))
if (toupper(ch@uplo) == "U") # which is TRUE
   crossprod(ch)
stopifnot(all.equal(as(crossprod(ch), "matrix"),
                    as(pm, "matrix"), tol=1e-14))
}
\keyword{classes}
\keyword{algebra}
