\name{schur}
\title{Schur Decomposition of a Matrix}
\usage{
schur.Matrix(x, vectors=TRUE, rcond=FALSE, threshold=0)
}
\alias{schur}
\alias{schur.Matrix}
\description{
    Computes the Schur decomposition and eigenvalues of a square matrix.
}
\arguments{
    \item{x}{
	numeric or complex square Matrix inheriting from class
	\code{"Matrix"}. Missing values (NAs) are not allowed.
    }
    \item{vectors}{logical.  When \code{TRUE} (the default), the Schur
	vectors are computed.
    }
    \item{rcond}{logical vector of length 2.  When \code{rcond[1]} is
	\code{TRUE}, the reciprocal condition numbers for the
	eigenvalues are computed.  When \code{rcond[2]} is
	\code{TRUE}, the reciprocal condition numbers for the
	eigenvalues are computed.  If a single logical value is given,
	its value is used for both conditions.  Default is \code{FALSE},
	indicating that neither set of reciprocal condition numbers is
	to be computed.
    }
    \item{threshold}{
	a numeric value. A diagonal block associated with a particular
	eigenvalue will appear in the top left portion of the Schur
	decomposition if it is greater than or equal to \code{threshold} in
	magnitude. By default, \code{threshold} is zero,  so that
	eigenvalues are unordered. Further ordering can be obtained by  
	applying \code{schurmod} to an existing Schur decomposition.
    }
%  \item{tune}{
%  a integer vector or list of named tuning parameters that may affect
%  computational efficiency. The relevant parameters are 
%  \code{NB}, \code{NBMIN}, \code{NX}, \code{NS}, and \code{MAXB},
% as described in \code{.laenv}.
%  The default is to use the settings in \code{.laenv}, which are initialized in
%  but may be changed by the user.
%  }
}
\value{
    An object of class \code{c("schur.Matrix", "decomp")} whose
    attributes include the eigenvalues and condition numbers (if requested).
}
\details{
    Based on the Lapack functions \code{dgeesx} and \code{zgeesx}.
%  Parameters in \code{tune} may affect performance and workspace requirements, 
%  but the optimal settings are both machine and problem dependent. 
%  Users may want to experiment with \code{tune} (see Anderson et al., p. 72-74) to 
%  obtain improved performance.
}
\section{BACKGROUND}{
    If \code{A} is a square matrix, then \code{A = Q T t(Q)}, where
    \code{Q} is orthogonal, and \code{T} is upper quasi-triangular
    (nearly triangular with either 1 by 1 or 2 by 2 blocks on the
    diagonal).
    The eigenvalues of \code{A} are the same as those of \code{T},
    which are easy to compute. The Schur form is used most often for
    computing non-symmetric eigenvalue decompositions, and for computing
    functions of matrices such as matrix exponentials.
}
\references{
Anderson, E., et al. (1994).
\emph{LAPACK User's Guide,}
2nd edition, SIAM, Philadelphia.
}
\seealso{
%    \code{\link{schur.Matrix.object}}, \code{\link{schurmod}},
%    \code{\link{eigen.Matrix}}, \code{\link{.laenv}}
}
\examples{
A <- Matrix(rnorm( 9*9, sd = 100), nrow = 9)
\dontrun{
schur.A <- schur(A)                          # schur factorization
mod.eig <- Mod(attr(schur.A, "eigenvalues")) # eigenvalue modulus
schur.A <- schurmod(schur.A)                 # reordered factorization
}
}
\keyword{algebra}
