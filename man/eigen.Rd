\name{eigen}
\title{Spectral Decomposition of a Matrix}
\usage{
eigen(x, symmetric, only.values = TRUE)
eigen.Matrix(x, vectors = TRUE, balance = FALSE, rcond = FALSE, schur = FALSE)
}
\alias{eigen}
\alias{eigen.default}
\alias{eigen.Matrix}
\arguments{
    \item{x}{a matrix whose spectral decomposition is to be computed.}
    \item{vectors}{if \code{TRUE}, both the eigenvalues and eigenvectors
	are returned, otherwise only the eigenvalues are
	returned.}
    \item{balance}{logical.  Not used at present.}
    \item{rcond}{logical.  Not used at present.}
    \item{schur}{logical.  Not used at present.}
}
\description{
    \code{eigen} is a generic function.  The default method
    \code{eigen.default} computes eigenvalues and eigenvectors by
    providing an interface to the EISPACK routines \code{RS}, \code{RG},
    \code{CH} and \code{CG}.  \code{eigen.Matrix} provides an interface
    to the Lapack functions \code{DSYEV} and \code{DGEEV}.
}
\value{
  The spectral decomposition of \code{x} is returned
  as components of a list.
  \item{values}{a vector containing the \eqn{p} eigenvalues of \code{x},
    sorted in \emph{decreasing} order, according to \code{Mod(values)}
    if they are complex.}
  \item{vectors}{a \eqn{p\times p}{p * p} matrix whose columns contain the
    eigenvectors of \code{x}, or \code{NULL} if \code{only.values} is
    \code{TRUE}.}
}
\references{
  Smith, B. T, Boyle, J. M., Dongarra, J. J., Garbow, B. S., Ikebe,Y.,
  Klema, V., and  Moler, C. B. (1976).
  \emph{Matrix Eigensystems Routines -- EISPACK Guide}.
  Springer-Verlag Lecture Notes in Computer Science.

  Anderson, E., et al. (1999). \emph{LAPACK User's Guide,}
  3rd edition, SIAM, Philadelphia.
}
\seealso{
  \code{\link{svd}}, a generalization of \code{eigen}; \code{\link{qr}}, and
  \code{\link{chol}} for related decompositions.

  To compute the determinant of a matrix, the \code{\link{qr}}
  decomposition is much more efficient: \code{\link{det}}.
}
\examples{
eigen(cbind(c(1,-1),c(-1,1)))   # uses Eispack
eigen(cbind(c(1,-1),c(-1,1)), symmetric = FALSE)# Eispack (different algorithm).
eigen(as.Matrix(cbind(c(1,-1),c(-1,1))))  # uses Lapack

eigen(cbind(1,c(1,-1)), only.values = TRUE)
eigen(as.Matrix(cbind(1,c(1,-1))), vectors = FALSE)
eigen(cbind(-1,2:1)) # complex values
eigen(as.Matrix(cbind(-1,2:1))) # complex values from Lapack
eigen(print(cbind(c(0,1i), c(-1i,0))))# Hermite ==> real Eigen values
## 3 x 3:
eigen(cbind( 1,3:1,1:3))
eigen(as.Matrix(cbind( 1,3:1,1:3)))
eigen(cbind(-1,c(1:2,0),0:2)) # complex values
eigen(as.Matrix(cbind(-1,c(1:2,0),0:2))) # complex values

Meps <- .Alias(.Machine$double.eps)
m <- matrix(round(rnorm(25),3), 5,5)
sm <- m + t(m) #- symmetric matrix
eigen(as.Matrix(sm), vectors = FALSE) # ordered INcreasingly
em <- eigen(sm); V <- em$vect
print(lam <- em$values) # ordered DEcreasingly


stopifnot(
 abs(sm \%*\% V - V \%*\% diag(lam))          < 60*Meps,
 abs(sm       - V \%*\% diag(lam) \%*\% t(V)) < 60*Meps)

##------- Symmetric = FALSE:  -- different to above : ---

em <- eigen(sm, symmetric = FALSE); V2 <- em$vect
print(lam2 <- em$values) # ordered decreasingly in ABSolute value !
                         # and V2 is not normalized (where V is):
print(i <- rev(order(lam2)))
stopifnot(abs(1 - lam2[i] / lam) < 60 * Meps)

zapsmall(Diag <- t(V2) \%*\% V2) # orthogonal, but not normalized
print(norm2V <- apply(V2 * V2, 2, sum))
stopifnot( abs(1- norm2V / diag(Diag)) < 60*Meps)

V2n <- sweep(V2,2, STATS= sqrt(norm2V), FUN="/")## V2n are now Normalized EV
apply(V2n * V2n, 2, sum)
##[1] 1 1 1 1 1

## Both are now TRUE:
stopifnot(abs(sm \%*\% V2n - V2n \%*\% diag(lam2))            < 60*Meps,
          abs(sm         - V2n \%*\% diag(lam2) \%*\% t(V2n)) < 60*Meps)

## Re-ordered as with symmetric:
sV <- V2n[,i]
slam <- lam2[i]
all(abs(sm \%*\% sV -  sV \%*\% diag(slam))             < 60*Meps)
all(abs(sm        -  sV \%*\% diag(slam) \%*\% t(sV)) < 60*Meps)
## sV  *is* now equal to V  -- up to sign (+-) and rounding errors
all(abs(c(1 - abs(sV / V)))       <     1000*Meps) # TRUE (P ~ 0.95)
}
\keyword{algebra}
\keyword{array}
