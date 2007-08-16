\name{Cholesky}
\docType{genericFunction}
\alias{Cholesky}
\alias{Cholesky,dsCMatrix-method}
\title{Cholesky Decomposition of a Sparse Matrix}
\usage{
Cholesky(A, perm = TRUE, LDL = TRUE, super = FALSE, \dots)
}
\description{
  Computes the Cholesky decomposition of a sparse, symmetric,
  positive-definite matrix.  However, typically \code{\link{chol}()}
  should rather be used unless you are interested in the different kinds
  of sparse Cholesky decompositions.
}
\arguments{
  \item{A}{sparse symmetric matrix.  No missing values or IEEE special
    values are allowed.}
  \item{perm}{logical scalar indicating if a fill-reducing permutation
    should be computed and applied to the rows and columns of \code{A}.
  Default is \code{TRUE}.}
  \item{LDL}{logical scalar indicating if the decomposition should be
    computed as LDL' where \code{L} is a unit lower triangular matrix.
    The alternative is LL' where \code{L} is lower triangular with
    arbitrary diagonal elements.  Default is \code{TRUE}.}
  \item{super}{logical scalar indicating is a supernodal decomposition
    should be created.  The alternative is a simplicial decomposition.
    Default is \code{FALSE}.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\value{
  an object inheriting from either
  \code{"\linkS4class{CHMsuper}"}, or
  \code{"\linkS4class{CHMsimpl}"}, depending on the \code{super}
  argument; both classes extend \code{"\linkS4class{CHMfactor}"} which
  extends \code{"\linkS4class{MatrixFactorization}"}.

  In other words, the result of \code{Cholesky()} is \emph{not} a
  matrix, and if you want one, you should probably rather use
  \code{\link{chol}()}.
}
\details{
  This is a generic function with special methods for different types
  of matrices.  Use \code{\link{showMethods}("Cholesky")} to list all
  the methods for the \code{\link{Cholesky}} generic.

  The method for class \code{\linkS4class{dsCMatrix}} of sparse matrices
  --- the only one available currently ---
  is based on functions from the CHOLMOD library.

  Again: If you just want the Cholesky decomposition of a matrix, you
  should probably rather use \code{\link{chol}(.)}.
}
\references{
  Tim Davis (2005)
  \emph{{CHOLMOD}: sparse supernodal {Cholesky} factorization and
    update/downdate}
  \url{http://www.cise.ufl.edu/research/sparse/cholmod/}

  Timothy A. Davis (2006)
  \emph{Direct Methods for Sparse Linear Systems}, SIAM Series
  \dQuote{Fundamentals of Algorithms}.
}
\seealso{
  Class definitions \code{\linkS4class{CHMfactor}} and
  \code{\linkS4class{dsCMatrix}} and function \code{\link{expand}}.

  Note that \code{\link{chol}()} returns matrices (inheriting from
  \code{"\linkS4class{Matrix}"}) whereas \code{Cholesky()} returns a
  \code{"\linkS4class{CHMfactor}"} object, and hence a typical user
  will rather use \code{chol(A)}.
}
\examples{
data(KNex)
mtm <- with(KNex, crossprod(mm))
str(mtm@factors) # empty list()
Cholesky(mtm)             # uses show(<MatrixFactorization>)
str(mtm@factors) # 'sPDCholesky' (simpl)
(Cm <- Cholesky(mtm, super = TRUE))
str(mtm@factors) # 'sPDCholesky'  *and*  'SPDCholesky'
str(cmat <- as(Cm, "sparseMatrix"))
cmat[1:20, 1:20]
b <- matrix(c(rep(0, 711), 1), nc = 1)
## solve(Cm, b) by default solves  Ax = b, where A = Cm'Cm !
stopifnot(identical(solve(Cm, b), solve(Cm, b, system = "A")),
          all.equal(solve(mtm, b),
                    solve(Cm, b)))
}
\keyword{array}
\keyword{algebra}
