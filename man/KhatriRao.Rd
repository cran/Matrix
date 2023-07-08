\name{KhatriRao}
\title{Khatri-Rao Matrix Product}
%
\keyword{algebra}
\keyword{arith}
\keyword{array}
\keyword{utilities}
%
\alias{KhatriRao}
%
\description{
  Computes Khatri-Rao products for any kind of matrices.

  The Khatri-Rao product is a column-wise Kronecker product.  Originally
  introduced by Khatri and Rao (1968), it has many different applications,
  see Liu and Trenkler (2008) for a survey.  Notably, it is used in
  higher-dimensional tensor decompositions, see Bader and Kolda (2008).
}
\usage{
KhatriRao(X, Y = X, FUN = "*", sparseY = TRUE, make.dimnames = FALSE)
}
\arguments{
  \item{X,Y}{matrices of with the same number of columns.}
  \item{FUN}{the (name of the) \code{\link{function}} to be used for the
    column-wise Kronecker products, see \code{\link{kronecker}},
    defaulting to the usual multiplication.}
  \item{sparseY}{logical specifying if \code{Y} should be coerced and
    treated as \code{\linkS4class{sparseMatrix}}.   Set this to
  \code{FALSE}, e.g., to distinguish structural zeros from zero entries.}
  \item{make.dimnames}{logical indicating if the result should inherit
    \code{\link{dimnames}} from \code{X} and \code{Y} in a simple way.}
}
%\details{}
\value{
  a \code{"\linkS4class{CsparseMatrix}"}, say \code{R}, the Khatri-Rao
  product of \code{X} (\eqn{n \times k}{n x k}) and \code{Y} (\eqn{m
    \times k}{m x k}), is of dimension \eqn{(n\cdot m) \times k}{(n*m) x k},
  where the j-th column, \code{R[,j]} is the kronecker product
  \code{\link{kronecker}(X[,j], Y[,j])}.
}
\note{%% TODO? Could make it generic, and have dense and sparse methods
  The current implementation is efficient for large sparse matrices.
}
\references{
  Khatri, C. G., and Rao, C. Radhakrishna (1968)
  Solutions to Some Functional Equations and Their Applications to
  Characterization of Probability Distributions.
  \emph{Sankhya: Indian J. Statistics, Series A} \bold{30}, 167--180.

  Bader, Brett W, and Tamara G Kolda (2008)
  Efficient MATLAB Computations with Sparse and Factored Tensors.
  \emph{SIAM J. Scientific Computing} \bold{30}, 205--231.
}
\author{
  Original by Michael Cysouw, Univ. Marburg;
  minor tweaks, bug fixes etc, by Martin Maechler.
}
\seealso{
 \code{\link{kronecker}}.
}
\examples{
## Example with very small matrices:
m <- matrix(1:12,3,4)
d <- diag(1:4)
KhatriRao(m,d)
KhatriRao(d,m)
dimnames(m) <- list(LETTERS[1:3], letters[1:4])
KhatriRao(m,d, make.dimnames=TRUE)
KhatriRao(d,m, make.dimnames=TRUE)
dimnames(d) <- list(NULL, paste0("D", 1:4))
KhatriRao(m,d, make.dimnames=TRUE)
KhatriRao(d,m, make.dimnames=TRUE)
dimnames(d) <- list(paste0("d", 10*1:4), paste0("D", 1:4))
(Kmd <- KhatriRao(m,d, make.dimnames=TRUE))
(Kdm <- KhatriRao(d,m, make.dimnames=TRUE))

nm <- as(m, "nsparseMatrix")
nd <- as(d, "nsparseMatrix")
KhatriRao(nm,nd, make.dimnames=TRUE)
KhatriRao(nd,nm, make.dimnames=TRUE)

stopifnot(dim(KhatriRao(m,d)) == c(nrow(m)*nrow(d), ncol(d)))
## border cases / checks:
zm <- nm; zm[] <- FALSE # all FALSE matrix
stopifnot(all(K1 <- KhatriRao(nd, zm) == 0), identical(dim(K1), c(12L, 4L)),
          all(K2 <- KhatriRao(zm, nd) == 0), identical(dim(K2), c(12L, 4L)))

d0 <- d; d0[] <- 0; m0 <- Matrix(d0[-1,])
stopifnot(all(K3 <- KhatriRao(d0, m) == 0), identical(dim(K3), dim(Kdm)),
	  all(K4 <- KhatriRao(m, d0) == 0), identical(dim(K4), dim(Kmd)),
	  all(KhatriRao(d0, d0) == 0), all(KhatriRao(m0, d0) == 0),
	  all(KhatriRao(d0, m0) == 0), all(KhatriRao(m0, m0) == 0),
	  identical(dimnames(KhatriRao(m, d0, make.dimnames=TRUE)), dimnames(Kmd)))

## a matrix with "structural" and non-structural zeros:
m01 <- new("dgCMatrix", i = c(0L, 2L, 0L, 1L), p = c(0L, 0L, 0L, 2L, 4L),
           Dim = 3:4, x = c(1, 0, 1, 0))
D4 <- Diagonal(4, x=1:4) # "as" d
DU <- Diagonal(4)# unit-diagonal: uplo="U"
(K5  <- KhatriRao( d, m01))
K5d  <- KhatriRao( d, m01, sparseY=FALSE)
K5Dd <- KhatriRao(D4, m01, sparseY=FALSE)
K5Ud <- KhatriRao(DU, m01, sparseY=FALSE)
(K6  <- KhatriRao(diag(3),     t(m01)))
K6D  <- KhatriRao(Diagonal(3), t(m01))
K6d  <- KhatriRao(diag(3),     t(m01), sparseY=FALSE)
K6Dd <- KhatriRao(Diagonal(3), t(m01), sparseY=FALSE)
stopifnot(exprs = {
    all(K5 == K5d)
    identical(cbind(c(7L, 10L), c(3L, 4L)),
              which(K5 != 0, arr.ind = TRUE, useNames=FALSE))
    identical(K5d, K5Dd)
    identical(K6, K6D)
    all(K6 == K6d)
    identical(cbind(3:4, 1L),
              which(K6 != 0, arr.ind = TRUE, useNames=FALSE))
    identical(K6d, K6Dd)
})
}
