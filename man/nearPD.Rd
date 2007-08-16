\name{nearPD}
\alias{nearPD}
\title{Nearest Matrix to a Positive Definite Matrix}
\description{
 Computes the nearest positive definite matrix to an approximate
 one, typically a correlation or variance-covariance matrix.
}
\usage{
nearPD(x, corr = FALSE,
       eig.tol = 1e-06, conv.tol = 1e-07, posd.tol = 1e-08,
       do2eigen = TRUE, maxit = 100, verbose = FALSE)
}
\arguments{
  \item{x}{numeric \eqn{n \times n}{n * n} approximately positive
    definite matrix, typically an approximation to a correlation or
    covariance matrix.}
  \item{corr}{logical indicating if the matrix should be a
    \emph{correlation} matrix.}
  \item{eig.tol}{defines relative positiveness of eigenvalues compared
    to largest one.}
  \item{conv.tol}{convergence tolerance for algorithm.}
  \item{posd.tol}{tolerance for enforcing positive definiteness.}
  \item{do2eigen}{logical indicating if a
    \code{\link[sfsmisc]{posdefify}()} eigen step should be applied to
    the result of the Hingham algorithm.}

  \item{maxit}{maximum number of iterations allowed.}
  \item{verbose}{logical; if \code{TRUE} the iterations are monitored by
    print out.}
}
\details{
  %% FIXME: add more

  Note that setting \code{corr = TRUE} just sets \code{diag(.) <- 1}
  within the algorithm.
}
\value{
  an S3 object of \code{\link{class}} \code{"nearPD"}, basically a list
  with components
  \item{mat}{a matrix of class \code{\linkS4class{dpoMatrix}}, the
  computed positive-definite matrix.}
  \item{corr}{logical, just the argument \code{corr}.}
  \item{normF}{the Frobenius norm (\code{\link{norm}(x-X, "F")}) of the
  difference between the original and the resulting matrix.}
  \item{iterations}{number of iterations needed.}
  \item{converged}{logical indicating if iterations converged.}
}
\references{%% more in /u/maechler/R/Pkgs/sfsmisc/man/posdefify.Rd
  Cheng, Sheung Hun and Higham, Nick (1998)
  A Modified Cholesky Algorithm Based on a Symmetric Indefinite Factorization;
  \emph{SIAM J. Matrix Anal.\ Appl.}, \bold{19}, 1097--1110.

  Highham (2002)
  Computing the nearest correlation matrix - a problem from finance;
  \emph{IMA Journal of Numerical Analysis} \bold{22}, 329--343.
}
\author{Jens Oehlschlaegel donated a first version.  Subsequent changes
  by the Matrix package authors.
}
\seealso{More simple versions with a similar purpose by
  \code{\link[sfsmisc]{posdefify}()}.
}
\examples{
 set.seed(27)
 m <- matrix(round(rnorm(25),2), 5, 5)
 m <- m + t(m)
 diag(m) <- pmax(0, diag(m)) + 1
 (m <- round(cov2cor(m), 2))

 str(near.m <- nearPD(m))
 round(near.m$mat, 2)
 norm(m - near.m$mat) # 1.102

 if(require("sfsmisc")) {
    m2 <- posdefify(m) # a simpler approach
    norm(m - m2)  # 1.185, i.e., slightly "less near"
 }
}
\keyword{algebra}
\keyword{array}
