\name{qr.Matrix}
\alias{qr.Matrix}
\title{
QR Decomposition of a Matrix
}
\description{
Computes the QR decomposition of a Matrix.
}
\usage{
qr.Matrix(x, pivot=F, tune=<see below>, workspace=<see below>)
}
\arguments{
\item{x}{
numeric or complex matrix. Missing values (NAs) are not allowed.
}
\item{pivot}{
a logical vector of length equal to the number of columns in \code{x}, indicating
whether or not each column is free for pivoting. The parameter \code{pivot} can 
also be a single logical value, indicating either that all or no columns 
are free for pivoting. By default, no column pivoting is done.
}
\item{tune}{
a integer vector or list of named tuning parameters that may affect 
computational efficiency when no pivoting is done. The relevant parameters
are the blocksize parameters \code{NB}, \code{NBMIN}, and \code{NX} as described in \code{.laenv}.
This parameter has no effect when pivoting is allowed.
The default is to use the settings in \code{.laenv}, which in S-Plus is initialized
to use an unblocked algorithm (\code{NB = 1}), but may be changed by the user.
}
\item{workspace}{
workspace provided to the underlying software.
The default is to use the optimum value according to \code{tune}.
The optimal workspace for the problem (for the given \code{tune} specifications) is
included as part of the output attributes.
}}
\value{
The QR decomposition of \code{x} is returned as
an object of class \code{c("qr.Matrix", "decomp")}.
A copy of the call to \code{qr.Matrix} is included as an attribute.
}
\details{
Based on the functions dgeqrf, dgeqpf, zgeqrf, and zgeqpf from LAPACK
(Anderson et al. 1994).
Parameters in \code{tune} may affect performance and workspace requirements, 
but the optimal settings are both machine and problem dependent. 
Users may want to experiment with \code{tune} (see Anderson et al., p. 72-74) 
to obtain improved performance.
}
\references{
Anderson, E., et al. (1994).
\emph{LAPACK User's Guide,}
2nd edition, SIAM, Philadelphia.
}
\seealso{
\code{\link{qr.Matrix.object}}, \code{\link{facmul.qr.Matrix}}, \code{\link{expand.qr.Matrix}}, \code{\link{.laenv}}
}
\examples{
 A <- Matrix( rnorm(60), nrow = 12, ncol = 5)
 qr(A, pivot = T)                           # QR factorization of A
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
