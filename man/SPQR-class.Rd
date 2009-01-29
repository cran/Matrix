\name{SPQR-class}
\docType{class}
\alias{SPQR-class}
\title{Class "SPQR" of SuiteSparse QR multifrontal sparse QR Decompositions}
\description{
  Objects of class \code{"SPQR"} are typically the result of
  \code{\link{spqr}}, one kind of sparse QR decompositions.

  Currently, \code{\link{qr}(<sparseMatrix>)} is more widely used and
  provides a different class, \code{\linkS4class{sparseQR}}.
}
\section{Slots}{
  \describe{
    \item{\code{Q}:}{Object of class \code{"\linkS4class{dgCMatrix}"},
      the \eqn{Q} matrix.}
    \item{\code{R}:}{Object of class \code{"dgCMatrix"}, the \eqn{R} matrix.}
    \item{\code{p}:}{integer vector providing the 0-index \bold{p}ermutation.}
    \item{\code{rank}:}{integer scalar providing the rank as
      \dQuote{revealed} by the algorithm.}
    \item{\code{Dim}:}{Object of class \code{"integer"}, the dimension
      of the original matrix.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{MatrixFactorization}"}, directly.
}
% \section{Methods}{
% No methods defined with class "SPQR" in the signature.
% }
% \references{ ~put references to the literature/web site here ~ }
\seealso{
  \code{\link{spqr}} which produces \code{"SPQR"} classed results.

  The (currently much wider used) alternative \code{\link{qr}} which
  produces \code{\linkS4class{sparseQR}} classed results.
}
\examples{
showClass("SPQR")
## See  the   example(spqr)  ones
}
\keyword{classes}
