\name{MatrixFactorization-class}
\title{Virtual Class "MatrixFactorization" of Matrix Factorizations}
%
\docType{class}
\keyword{algebra}
\keyword{array}
\keyword{classes}
%
\alias{MatrixFactorization-class}
%
\alias{BunchKaufmanFactorization-class}
\alias{CholeskyFactorization-class}
\alias{SchurFactorization-class}
\alias{LU-class}
\alias{QR-class}
%
\alias{determinant,MatrixFactorization,missing-method}
\alias{dim,MatrixFactorization-method}
\alias{dimnames,MatrixFactorization-method}
\alias{dimnames<-,MatrixFactorization,NULL-method}
\alias{dimnames<-,MatrixFactorization,list-method}
\alias{length,MatrixFactorization-method}
\alias{show,MatrixFactorization-method}
\alias{unname,MatrixFactorization-method}
%
\alias{show,BunchKaufmanFactorization-method}
\alias{show,CholeskyFactorization-method}
\alias{show,SchurFactorization-method}
\alias{show,LU-method}
\alias{show,QR-method}
%
\description{
  \code{MatrixFactorization} is the virtual class of
  factorizations of \eqn{m \times n}{m-by-n} matrices \eqn{A},
  having the general form
  \deqn{P_{1} A P_{2} = A_{1} \cdots A_{p}}{P1 * A * P2 = A1 * ... * Ap}
  or (equivalently)
  \deqn{A = P_{1}' A_{1} \cdots A_{p} P_{2}'}{A = P1' * A1 * ... * Ap * P2'}
  where \eqn{P_{1}}{P1} and \eqn{P_{2}}{P2} are permutation matrices.
  Factorizations requiring symmetric \eqn{A} have the constraint
  \eqn{P_{2} = P_{1}'}{P2 = P1'}, and factorizations without row
  or column pivoting have the constraints
  \eqn{P_{1} = I_{m}}{P1 = Im} and \eqn{P_{2} = I_{n}}{P2 = In},
  where \eqn{I_{m}}{Im} and \eqn{I_{n}}{In} are the
  \eqn{m \times m}{m-by-m} and \eqn{n \times n}{n-by-n} identity matrices.
  
  \code{CholeskyFactorization}, \code{BunchKaufmanFactorization},
  \code{SchurFactorization}, \code{LU}, and \code{QR} are the virtual
  subclasses of \code{MatrixFactorization} containing all Cholesky,
  Bunch-Kaufman, Schur, LU, and QR factorizations, respectively.
}
\section{Slots}{
  \describe{
    \item{\code{Dim}}{an integer vector of length 2 giving the
      dimensions of the factorized matrix.}
    \item{\code{Dimnames}}{a list of length 2 preserving the
      \code{dimnames} of the factorized matrix.  Each element
      must be \code{NULL} or a character vector of length equal
      to the corresponding element of \code{Dim}.}
  }
}
\section{Methods}{
  \describe{
    \item{\code{determinant}}{\code{signature(x = "MatrixFactorization", logarithm = "missing")}:
      sets \code{logarithm = TRUE} and recalls the generic function.}
    \item{\code{dim}}{\code{signature(x = "MatrixFactorization")}:
      returns \code{x@Dim}.}
    \item{\code{dimnames}}{\code{signature(x = "MatrixFactorization")}:
      returns \code{x@Dimnames}.}
    \item{\code{dimnames<-}}{\code{signature(x = "MatrixFactorization", value = "NULL")}:
      returns \code{x} with \code{x@Dimnames} set to \code{list(NULL, NULL)}.}
    \item{\code{dimnames<-}}{\code{signature(x = "MatrixFactorization", value = "list")}:
      returns \code{x} with \code{x@Dimnames} set to \code{value}.}
    \item{\code{length}}{\code{signature(x = "MatrixFactorization")}:
      returns \code{prod(x@Dim)}.}
    \item{\code{show}}{\code{signature(object = "MatrixFactorization")}:
      prints the internal representation of the factorization using
      \code{\link{str}}.}
    \item{\code{solve}}{\code{signature(a = "MatrixFactorization", b = .)}:
      see \code{\link{solve-methods}}.}
    \item{\code{unname}}{\code{signature(obj = "MatrixFactorization")}:
      returns \code{obj} with \code{obj@Dimnames} set to
      \code{list(NULL, NULL)}.}
  }
}
\seealso{
  The virtual class \code{\linkS4class{compMatrix}}
  of \emph{factorizable} matrices.
  
  Classes extending \code{CholeskyFactorization}, namely
  \code{\linkS4class{Cholesky}}, \code{\linkS4class{pCholesky}},
  and \code{\linkS4class{CHMfactor}}.

  Classes extending \code{BunchKaufmanFactorization}, namely
  \code{\linkS4class{BunchKaufman}} and \code{\linkS4class{pBunchKaufman}}.

  Classes extending \code{SchurFactorization}, namely
  \code{\linkS4class{Schur}}.

  Classes extending \code{LU}, namely
  \code{\linkS4class{denseLU}} and \code{\linkS4class{sparseLU}}.
  
  Classes extending \code{QR}, namely \code{\linkS4class{sparseQR}}.
  
  Generic functions \code{\link{Cholesky}}, \code{\link{BunchKaufman}},
  \code{\link{Schur}}, \code{\link{lu}}, and \code{\link{qr}} for
  \emph{computing} factorizations.

  Generic functions \code{\link{expand1}} and \code{\link{expand2}}
  for constructing matrix factors from \code{MatrixFactorization}
  objects.
}
\examples{
showClass("MatrixFactorization")
}
