\name{Cholesky}
\docType{genericFunction}
\alias{Cholesky}
\alias{Cholesky,dsCMatrix-method}
\title{Cholesky Decomposition of a Sparse Matrix}
\usage{
Cholesky(A, perm, LDL, super, \dots)
}
\description{
  Computes the Cholesky decomposition of a sparse, symmetric,
  positive-definite matrix.
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
  an object of class \code{"\linkS4class{CHMfactor}"}.
}
\details{
  This is a generic function with special methods for different types
  of matrices.  Use \code{\link{showMethods}("Cholesky")} to list all
  the methods for the \code{\link{Cholesky}} generic.

  The method for class \code{\linkS4class{dsCMatrix}} of sparse matrices
  is based on functions from the CHOLMOD library.
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
  Class definitions \code{\linkS4class{CHMfactor}} and \code{\linkS4class{dsCMatrix}}
  and function \code{\link{expand}}.

  Note that \code{\link{chol}()} returns matrices (inheriting from
  \code{"\linkS4class{Matrix}"}) whereas \code{Cholesky()} returns a
  \code{"\linkS4class{CHMfactor}"} object.
}
\examples{
data(KNex)
mtm <- with(KNex, crossprod(mm))
Cholesky(mtm)             # uses show(<MatrixFactorization>)
Cholesky(mtm, super = TRUE)
}
\keyword{array}
\keyword{algebra}
