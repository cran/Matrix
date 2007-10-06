\name{lu}
\title{Triangular Decomposition of a Square Matrix}
\usage{
lu(x, \dots)
}
\alias{lu}
\description{
  Computes triangular decompositions of square matrices.
}
\arguments{
  \item{x}{a matrix. No missing values or IEEE special values are allowed.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\value{
  an object of class \code{"LU"}, i.e., \code{"\linkS4class{denseLU}"} or
  \code{"sparseLU"}, see \code{\linkS4class{sparseLU}}; this is
  a representation of a triangular decomposition of \code{x}.
}
\details{
  This is a generic function with special methods for different types
  of matrices.  Use \code{\link{showMethods}("lu")} to list all the methods
  for the \code{\link{lu}} generic.

  The method for class \code{\linkS4class{dgeMatrix}} is based on
  LAPACK's \code{"dgetrf"} subroutine.

  The method for class \code{\linkS4class{dgCMatrix}} of sparse matrices
  is based on functions from the CSparse library.
}
\references{
  Golub, G., and Van Loan, C. F. (1989).
  \emph{Matrix Computations,}
  2nd edition, Johns Hopkins, Baltimore.

  Tim Davis (2005)
  \url{http://www.cise.ufl.edu/research/sparse/CSparse/}

  Timothy A. Davis (2006)
  \emph{Direct Methods for Sparse Linear Systems}, SIAM Series
  \dQuote{Fundamentals of Algorithms}.
}
\seealso{
  Class definitions \code{\linkS4class{LU}} and \code{\linkS4class{sparseLU}}
  and function \code{\link{expand}};
  \code{\link{qr}}, \code{\link{chol}}.
}
\examples{

##--- Dense  -------------------------
x <- Matrix(rnorm(9), 3, 3)
lu(x)

##--- Sparse ------------------------

pm <- as(readMM(system.file("external/pores_1.mtx",
                            package = "Matrix")),
         "CsparseMatrix")
str(pmLU <- lu(pm))		# p is a 0-based permutation of the rows
                                # q is a 0-based permutation of the columns
## permute rows and columns of original matrix
ppm <- pm[pmLU@p + 1L, pmLU@q + 1L]
pLU <- pmLU@L \%*\% pmLU@U
## equal up to "rounding"
ppm[1:14, 1:5]
pLU[1:14, 1:5]  # product can have extra zeros
## "prove" consistency (up to rounding):
i0 <- ppm != pLU & ppm == 0
iN <- ppm != pLU & ppm != 0
stopifnot(all(abs((ppm - pLU)[i0]) < 1e-7), # absolute error for true 0
          all(abs((ppm - pLU)[iN]/ppm[iN]) < 1e-9)) # relative error
}
\keyword{array}
\keyword{algebra}
