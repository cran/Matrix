\name{SparseM-conversions}
%% only as long as not requiring R >= 3.2.0:
\newcommand{\CRANpkg}{\href{https://CRAN.R-project.org/package=#1}{\pkg{#1}}}
\title{Sparse Matrix Coercion from and to those from package \pkg{SparseM}}
\docType{methods}
\alias{SparseM-coerce-methods}
%
\alias{coerce,matrix.coo,CsparseMatrix-method}
\alias{coerce,matrix.coo,TsparseMatrix-method}
\alias{coerce,matrix.csc,CsparseMatrix-method}
\alias{coerce,matrix.csc,TsparseMatrix-method}
\alias{coerce,matrix.csr,CsparseMatrix-method}
\alias{coerce,matrix.csr,RsparseMatrix-method}
\alias{coerce,matrix.csr,TsparseMatrix-method}
%
\alias{coerce,CsparseMatrix,matrix.csr-method}
\alias{coerce,CsparseMatrix,matrix.coo-method}
\alias{coerce,CsparseMatrix,matrix.csc-method}
\alias{coerce,dsparseMatrix,matrix.csr-method}
\alias{coerce,dgRMatrix,matrix.csr-method}
\alias{coerce,dgCMatrix,matrix.csc-method}
\alias{coerce,dgTMatrix,matrix.coo-method}
\alias{coerce,matrix.csr,dgRMatrix-method}
\alias{coerce,matrix.csc,dgCMatrix-method}
\alias{coerce,matrix.coo,dgTMatrix-method}
\alias{coerce,matrix.csr,dgCMatrix-method}
\alias{coerce,matrix.coo,dgCMatrix-method}
%
\alias{coerce,matrix.csr,Matrix-method}
\alias{coerce,matrix.coo,Matrix-method}
\alias{coerce,matrix.csc,Matrix-method}
%
\description{
  Methods for coercion from and to sparse matrices from package \pkg{SparseM}
  are provided here, for ease of porting functionality to the
  \pkg{Matrix} package, and comparing functionality of the two
  packages.  All these work via the usual \code{\link{as}(., "<class>")}
  coercion, \preformatted{  as(from, Class)
}%pre
}
\section{Methods}{
  \describe{
    \item{from = "matrix.csr", to = "dgRMatrix"}{ ... }
    \item{from = "matrix.csc", to = "dgCMatrix"}{ ... }
    \item{from = "matrix.coo", to = "dgTMatrix"}{ ... }

    \item{from = "dgRMatrix", to = "matrix.csr"}{ ... }
    \item{from = "dgCMatrix", to = "matrix.csc"}{ ... }
    \item{from = "dgTMatrix", to = "matrix.coo"}{ ... }

    \item{from = "sparseMatrix", to = "matrix.csr"}{ ... }
    \item{from = "matrix.csr", to = "dgCMatrix"}{ ... }
    \item{from = "matrix.coo", to = "dgCMatrix"}{ ... }

    \item{from = "matrix.csr", to = "Matrix"}{ ... }
    \item{from = "matrix.csc", to = "Matrix"}{ ... }
    \item{from = "matrix.coo", to = "Matrix"}{ ... }
  }
}
\seealso{
  The documentation in CRAN package \CRANpkg{SparseM}, such as
  \code{\link[SparseM]{SparseM.ontology}}, and one important class,
  \code{\link[SparseM:matrix.csr-class]{matrix.csr}}.
}
\keyword{methods}
