\name{qr-methods}
\docType{methods}
\alias{qr-methods}
\alias{qr,dgCMatrix-method}
\title{Methods for the qr generic} 
\description{
  The \code{"Matrix"} package provides methods for the QR decomposition
  of  special classes of matrices.  The generic function is
  \code{\link[base]{qr}} but methods defined in this package can take
  extra arguments.  In particular there is an option for determining a
  fill-reducing permutation of the columns of a sparse, rectangular matrix.
}
\section{Methods}{
  \describe{
    \item{x = "dgCMatrix"}{QR decomposition of a general sparse
      double-precision matrix with \code{nrow(x) >= ncol(x)}.  Returns
      an object of class \code{"sparseQR"}, see
      \code{\linkS4class{sparseQR}}.
    }
  }
}
\seealso{\code{\link[base]{qr}}, \code{\linkS4class{sparseQR}},
    \code{\linkS4class{dgCMatrix}}}
\keyword{methods}
\keyword{algebra}
\keyword{array}
