\name{expand}
\alias{expand}
\title{
Expand a Decomposition into Factors --- Generic Function
}
\description{
Expands decompositions stored in compact form into factors.
This is a generic function.
Functions with names beginning in \code{expand.} will be
methods for this function.
\code{"lu.Matrix"}, \code{"lu.Hermitian"}, \code{"qr.Matrix"}, \code{"eigen.Matrix"}.
}
\usage{
expand(x, \dots)
}
\arguments{
\item{x}{
a matrix decomposition. 
}
\item{...}{
the method for \code{"qr.Matrix"} has an additional argument.
}}
\value{
the expanded decomposition.
}
\section{NOTE}{
Factors for decompositions such as \code{lu} and \code{qr} can be stored in a compact
form. The function \code{expand} allows all factors to be fully expanded. 
}
\references{
Golub, G., and Van Loan, C. F. (1989).
\emph{Matrix Computations,}
2nd edition, Johns Hopkins, Baltimore.
}
\seealso{
\code{\link{expand.lu.Matrix}}, \code{\link{expand.lu.Hermitian}}, \code{\link{expand.qr.Matrix}}, 
\code{\link{expand.eigen.Matrix}}, \code{\link{facmul}}.
}
\examples{
library(Matrix)
x <- Matrix( 1:9, 3, 3)
expand(qr(x))
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
