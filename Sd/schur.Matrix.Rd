\name{schur.Matrix}
\alias{schur.Matrix}
\title{
Schur Decomposition of a Matrix
}
\description{
Computes the Schur decomposition and eigenvalues of a square matrix.
}
\usage{
schur.Matrix(x, vectors=T, rcond=F, threshold=0, tune=<see below>, 
             workspace=<see below>)
}
\arguments{
\item{x}{
numeric or complex square Matrix inheriting from class \code{"Matrix"}. Missing values (NAs) are not allowed.
}
\item{vectors}{
a logical value telling whether or not to compute the Schur vectors.
The default is to compute them.
}
\item{rcond}{
either a logical value or a logical vector of length 2 telling whether or 
not to compute reciprocal condition numbers for the eigenvalues and Schur
vectors, respectively, and return them as attributes. 
By default, neither is computed.
}
\item{threshold}{
a numeric value. A diagonal block associated with a particular eigenvalue will
appear in the top left portion of the Schur decomposition if it is greater 
than or equal to \code{threshold} in magnitude. By default, \code{threshold} is zero, 
so that eigenvalues are unordered. Further ordering can be obtained by 
applying \code{schurmod} to an existing Schur decomposition.
}
\item{tune}{
a integer vector or list of named tuning parameters that may affect
computational efficiency. The relevant parameters are 
\code{NB}, \code{NBMIN}, \code{NX}, \code{NS}, and \code{MAXB}, as described in \code{.laenv}.
The default is to use the settings in \code{.laenv}, which are initialized in
S-Plus but may be changed by the user.
}
\item{workspace}{
workspace provided to the underlying software.
The default is to use the minimum allowed value.
The optimal workspace for the problem (for the given \code{tune} specifications) is
included as part of the output attributes.
}}
\value{
The result is an object of class \code{c("schur.Matrix", "decomp")}.
Attributes include eigenvalues, condition numbers (if requested via \code{rcond}), 
and workspace information, and a copy of the call to \code{schur.Matrix}.
}
\details{
Based on the functions dgeesx and zgeesx from LAPACK (Anderson et al. 1994).
Parameters in \code{tune} may affect performance and workspace requirements, 
but the optimal settings are both machine and problem dependent. 
Users may want to experiment with \code{tune} (see Anderson et al., p. 72-74) to 
obtain improved performance.
}
\section{BUGS}{
The condition numbers are not computed correctly. A bug report has been filed
with LAPACK.
}
\references{
Anderson, E., et al. (1994).
\emph{LAPACK User's Guide,}
2nd edition, SIAM, Philadelphia.
}
\seealso{
\code{\link{schur.Matrix.object}}, \code{\link{schurmod}}, \code{\link{eigen.Matrix}}, \code{\link{.laenv}}
}
\examples{
 n <- 9
 A <- Matrix( rnorm(n*n,sd = 100), nrow = n, ncol = n)
 schur.A <- schur(A)                                 # schur factorization
 mod.eig <- Mod(attr( schur.A, "eigenvalues"))       # eigenvalue modulus
 schur.A <- schur(schur.A)                           # reordered factorization
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
