\name{svd.Matrix}
\alias{svd.Matrix}
\title{
    Singular Value Decomposition of a Matrix
}
\description{
    Computes singular values and (optionally) singular vectors of a matrix.
}
\usage{
svd.Matrix(x, vectors, tune, workspace)
}
\arguments{
    \item{x}{
	numeric or complex Matrix inheriting from class
	\code{"Matrix"}. NAs not allowed.
    }
    \item{vectors}{
	either a logical value telling whether or not to compute
	singular vectors, a logical vector of length 2 telling whether
	or not to compute the left and right singular vectors,
	respectively, or a numeric vector indicating the number of
	singular vectors to be computed.
	If this number is strictly between 0 and the corresponding
	dimension of \code{x}, then \code{min(dim(x))} singular vectors
	will be computed.
	The default is to return \code{min(dim(x))} left singular
	vectors and \code{min(dim(x))} right singular vectors.
    }
    \item{tune}{
	a integer vector or list of named tuning parameters that may
	affect computational efficiency. The relevant parameters are
	\code{NB}, \code{NBMIN}, \code{NX}, and \code{NXSVD}, as
	described in \code{.laenv}. 
	The default is to use the settings in \code{.laenv}, which are
	initialized in S-Plus but may be changed by the user.
    }
    \item{workspace}{
	workspace provided to the underlying software. 
	The default is to use the minimum allowed value.
	The optimal workspace for the problem (for the given \code{tune}
	specifications) is included as part of the output attributes.
    }
}
\value{
    The result is an object of class \code{c( "svd.Matrix", "decomp")}.
    Attributes include the dimensions of the underlying matrix
    and a flag indicating whether or not the matrix is complex
    (this information is lost if no singular vectors are requested),
    as well as workspace information and a copy of the call to \code{svd.Matrix}.
}
\details{
    Based on the functions dgesvd and zgesvd from LAPACK
    (Anderson et al. 1994).
    In order for \code{x} to be recovered from the SVD, the number of
    left and right singular vectors need not exceed the minimum
    dimension of \code{x}.
    Parameters in \code{tune} may affect performance and workspace
    requirements, but the optimal settings are both machine and problem
    dependent.
    Users may want to experiment with \code{tune} (see Anderson et al.,
    p. 72-74) to obtain improved performance.
}
\section{BUGS}{
    The alleged optimal workspace can be smaller than the minimum allowed
    workspace.
}
\references{
    Anderson, E., et al. (1994).
    \emph{LAPACK User's Guide,}
    2nd edition, SIAM, Philadelphia.
}
\seealso{
    \code{\link{svd.Matrix.object}}, \code{\link{qr.Matrix}},
    \code{\link{eigen.Matrix}}
}
\examples{
 A <- Matrix( rnorm(60), nrow = 12, ncol = 5)
 svd.A <- svd(A)                        # SVD of A
 rcond(svd.A)                           # reciprocal condition number of A
 b <- rnorm(12)
 solve(svd.A,b)                         # least-squares solution to Ax = b
}
\keyword{algebra}
% Converted by Sd2Rd version 1.10.
