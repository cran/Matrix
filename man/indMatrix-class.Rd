\name{indMatrix-class}
\title{Index Matrices}
%
\docType{class}
\keyword{array}
\keyword{classes}
%
\alias{indMatrix-class}
%
\alias{!,indMatrix-method}
\alias{-,indMatrix,missing-method}
\alias{Math,indMatrix-method}
\alias{Summary,indMatrix-method}
\alias{coerce,indMatrix,pMatrix-method}
\alias{coerce,list,indMatrix-method}
\alias{coerce,matrix,indMatrix-method}
\alias{coerce,numeric,indMatrix-method}
\alias{determinant,indMatrix,logical-method}
\alias{diag,indMatrix-method}
\alias{diag<-,indMatrix-method}
\alias{log,indMatrix-method}
\alias{t,indMatrix-method}
\alias{which,indMatrix-method}
%
\description{
  The \code{indMatrix} class is the class of row and column
  \emph{index} matrices, stored as 1-based integer index vectors.
  A row (column) index matrix is a matrix whose rows (columns)
  are standard unit vectors.  Such matrices are useful
  when mapping observations to discrete sets of covariate values.

  Multiplying a matrix on the left by a row index matrix is
  equivalent to indexing its rows, i.e., sampling the rows
  \dQuote{with replacement}.  Analogously, multiplying a matrix
  on the right by a column index matrix is equivalent to
  indexing its columns.  Indeed, such products are implemented in
  \pkg{Matrix} as indexing operations; see \sQuote{Details} below.

  A matrix whose rows \emph{and} columns are standard unit vectors
  is called a \emph{permutation} matrix.  This special case is
  designated by the \code{\linkS4class{pMatrix}} class, a direct
  subclass of \code{indMatrix}.
}
\section{Objects from the Class}{
  Objects can be created explicitly with calls of the form
  \code{new("indMatrix", ...)}, but they are more commonly created
  by coercing 1-based integer index vectors, with calls of the
  form \code{as(., "indMatrix")}; see \sQuote{Methods} below.
}
\section{Slots}{
  \describe{
    \item{\code{margin}}{an integer, either 1 or 2, specifying
      whether the matrix is a row (1) or column (2) index.}
    \item{\code{perm}}{a 1-based integer index vector, i.e.,
      a vector of length \code{Dim[margin]} with elements
      taken from \code{1:Dim[1+margin\%\%2]}.}
    \item{\code{Dim},\code{Dimnames}}{inherited from virtual
      superclass \code{\linkS4class{Matrix}}.}
  }
}
\section{Extends}{
  Classes \code{"\linkS4class{sparseMatrix}"} and
  \code{"\linkS4class{generalMatrix}"}, directly.
}
\section{Methods}{
  \describe{
    \item{\code{\%*\%}}{\code{signature(x = "indMatrix", y = "Matrix")}
      and others listed by \code{showMethods("\%*\%", classes = "indMatrix")}:
      matrix products implemented where appropriate as indexing operations.}
    \item{\code{coerce}}{\code{signature(from = "numeric", to = "indMatrix")}:
      supporting typical \code{indMatrix} construction from
      a vector of positive integers.  Row indexing is assumed.}
    \item{\code{coerce}}{\code{signature(from = "list", to = "indMatrix")}:
      supporting \code{indMatrix} construction for row \emph{and}
      column indexing, including index vectors of length 0 and
      index vectors whose maximum is less than the number of rows
      or columns being indexed.}
    \item{\code{coerce}}{\code{signature(from = "indMatrix", to = "matrix")}:
      coercion to a traditional \code{\link{matrix}} of \link{logical} type,
      with \code{FALSE} and \code{TRUE} in place of 0 and 1.}
    \item{\code{t}}{\code{signature(x = "indMatrix")}:
      the transpose, which is an \code{indMatrix} with identical
      \code{perm} but opposite \code{margin}.}
    \item{\code{rowSums},\code{rowMeans},\code{colSums},\code{colMeans}}{\code{signature(x = "indMatrix")}:
      row and column sums and means.}
    \item{\code{rbind2},\code{cbind2}}{\code{signature(x = "indMatrix", y = "indMatrix")}:
      row-wise catenation of two row index matrices with equal numbers
      of columns and column-wise catenation of two column index matrices
      with equal numbers of rows.}
    \item{kronecker}{\code{signature(X = "indMatrix", Y = "indMatrix")}:
      Kronecker product of two row index matrices or two column index
      matrices, giving the row or column index matrix corresponding to
      their \dQuote{interaction}.}
  }
}
\author{Fabian Scheipl and Uni Muenchen, building on the existing class
  \code{\linkS4class{pMatrix}} after a nice hike's conversation with
  Martin Maechler.  Methods for \code{\link{crossprod}(x, y)} and
  \code{\link{kronecker}(x, y)} with both arguments inheriting from
  \code{indMatrix} were made considerably faster thanks to a suggestion
  by Boris Vaillant.  Diverse tweaks by Martin Maechler and
  Mikael Jagan, notably the latter's implementation of \code{margin},
  prior to which the \code{indMatrix} class was designated only for
  row index matrices.
}
\details{
  The transpose of an index matrix is an index matrix with identical
  \code{perm} but opposite \code{margin}.  Hence the transpose of a
  row index matrix is a column index matrix, and vice versa.

  The cross product of a row index matrix \code{R} and itself is a
  diagonal matrix whose diagonal entries are the the number of entries
  in each column of \code{R}.
  
  Given a row index matrix \code{R} with \code{perm} slot \code{p},
  a column index matrix \code{C} with \code{perm} slot \code{q},
  and a matrix \code{M} with conformable dimensions, we have
  \tabular{lclcl}{
    \eqn{R M} \tab = \tab \code{R \%*\% M}        \tab = \tab \code{M[p, ]}\cr
    \eqn{M C} \tab = \tab \code{M \%*\% C}        \tab = \tab \code{M[, q]}\cr
    \eqn{C'M} \tab = \tab \code{crossprod(C, M)}  \tab = \tab \code{M[q, ]}\cr
    \eqn{MR'} \tab = \tab \code{tcrossprod(M, R)} \tab = \tab \code{M[, p]}\cr
    \eqn{R'R} \tab = \tab \code{crossprod(R)}     \tab = \tab \code{Diagonal(x=tabulate(p, ncol(R)))}\cr
    \eqn{CC'} \tab = \tab \code{tcrossprod(C)}    \tab = \tab \code{Diagonal(x=tabulate(q, nrow(C)))}
  }

  Operations on index matrices that result in index matrices will
  accordingly return an \code{indMatrix}.  These include products
  of two column index matrices and (equivalently) column-indexing
  of a column index matrix (when dimensions are not dropped).
  Most other operations on \code{indMatrix} treat them as sparse
  nonzero pattern matrices (i.e., inheriting from virtual class
  \code{\linkS4class{nsparseMatrix}}).  Hence vector-valued subsets
  of \code{indMatrix}, such as those given by \code{\link{diag}},
  are always of type \code{"\link{logical}"}.
}
\seealso{
  Subclass \code{\linkS4class{pMatrix}} of permutation matrices,
  a special case of index matrices; virtual class
  \code{\linkS4class{nMatrix}} of nonzero pattern matrices,
  and its subclasses.
}
\examples{
p1 <- as(c(2,3,1), "pMatrix")
(sm1 <- as(rep(c(2,3,1), e=3), "indMatrix"))
stopifnot(all(sm1 == p1[rep(1:3, each=3),]))

## row-indexing of a <pMatrix> turns it into an <indMatrix>:
class(p1[rep(1:3, each=3),])

set.seed(12) # so we know '10' is in sample
## random index matrix for 30 observations and 10 unique values:
(s10 <- as(sample(10, 30, replace=TRUE),"indMatrix"))

## Sample rows of a numeric matrix :
(mm <- matrix(1:10, nrow=10, ncol=3))
s10 \%*\% mm

set.seed(27)
IM1 <- as(sample(1:20, 100, replace=TRUE), "indMatrix")
IM2 <- as(sample(1:18, 100, replace=TRUE), "indMatrix")
(c12 <- crossprod(IM1,IM2))
## same as cross-tabulation of the two index vectors:
stopifnot(all(c12 - unclass(table(IM1@perm, IM2@perm)) == 0))

# 3 observations, 4 implied values, first does not occur in sample:
as(2:4, "indMatrix")
# 3 observations, 5 values, first and last do not occur in sample:
as(list(2:4, 5), "indMatrix")

as(sm1, "nMatrix")
s10[1:7, 1:4] # gives an "ngTMatrix" (most economic!)
s10[1:4, ]  # preserves "indMatrix"-class

I1 <- as(c(5:1,6:4,7:3), "indMatrix")
I2 <- as(7:1, "pMatrix")
(I12 <- rbind(I1, I2))
stopifnot(is(I12, "indMatrix"),
          identical(I12, rbind(I1, I2)),
	  colSums(I12) == c(2L,2:4,4:2))
}
