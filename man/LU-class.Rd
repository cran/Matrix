\name{LU-class}
\docType{class}
\alias{LU-class}
\alias{expand,LU-method}
\title{LU Matrix Decompositions}
\description{
  The \code{"LU"} class is the class of LU decompositions of
  real matrices.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("LU", ...)}.
  More commonly the objects are created explicitly from calls of the form
  \code{lu(mm)} where \code{mm} is an object that inherits from the
  \code{"dgeMatrix"} class or as a side-effect of other functions applied
  to \code{"dgeMatrix"} objects.
}
\section{Slots}{
  \describe{
    \item{\code{x}:}{object of class \code{"numeric"}.  The \code{"L"}
      (unit lower triangular) and \code{"U"} (upper triangular) factors
      of the original matrix.  These are stored in a packed format
      described in the Lapack manual.}
    \item{\code{Dim}:}{Object of class \code{"integer"} - the dimensions
     of the original matrix - must be an integer vector with exactly two
     non-negative values.}
    \item{\code{pivot}:}{Object of class \code{"integer"} - a vector of
      length \code{min(Dim)} that describes the permutation applied to
      the rows of the original matrix.  The contents of this vector are
      described in the Lapack manual.}
  }
}
\section{Methods}{
  \describe{
    \item{expand}{\code{signature(x = "LU")}: Produce the \code{"L"} and
      \code{"U"} factors as a named list of matrices.}
  }
}
% \references{}
% \author{}
% \note{}
\seealso{
  \code{\link{dgeMatrix-class}}, \code{\link{lu}}, \code{\link{expand}}
}
\examples{
mm <- Matrix(round(rnorm(9),2), nrow = 3, ncol = 3)
mm
str(lum <- lu(mm))
elu <- expand(lum)
elu
elu$L \%*\% elu$U
}
\keyword{classes}
\keyword{algebra}
