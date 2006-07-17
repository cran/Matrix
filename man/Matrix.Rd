\name{Matrix}
\alias{Matrix}
\title{Construct a Classed Matrix}
\usage{
Matrix(data=NA, nrow=1, ncol=1, byrow=FALSE, dimnames=NULL,
       sparse = NULL)
}
\description{
  Construct a Matrix of a class that inherits from \code{Matrix}.
}
\arguments{
  \item{data}{an optional numeric data vector or matrix.}
  \item{nrow}{the desired number of rows}
  \item{ncol}{the desired number of columns}
  \item{byrow}{logical.  If \code{FALSE} (the default) the matrix is
    filled by columns, otherwise the matrix is filled by rows.}
  \item{dimnames}{a \code{\link{dimnames}} attribute for the matrix: a
    \code{list} of two character components.  They are set if not
    \code{\link{NULL}} (as per default).}
  \item{sparse}{logical or \code{NULL}, specifying if the result should
    be sparse or not.  By default, it is made sparse when more than half
    of the entries are 0.}
}
\value{
  Returns an \code{nrow} by \code{ncol} matrix of a class that
  inherits from \code{"Matrix"}.
}
\details{
  If either of \code{nrow} or \code{ncol} is not given, an attempt is
  made to infer it from the length of \code{data} and the other
  parameter.
  Further, \code{Matrix()} makes efforts to keep \code{\link{logical}}
  matrices logical, i.e., inheriting from class \code{\linkS4class{lMatrix}},
  and to determine specially structured matrices such as symmetric,
  triangular or diagonal ones.  Note that a \emph{symmetric} matrix also
  needs symmetric \code{\link{dimnames}}, e.g., by specifying
  \code{dimnames = list(NULL,NULL)}, see the examples.

  Most of the time, the function works via a traditional (\emph{full})
  \code{\link{matrix}}.  However, \code{Matrix(0, nrow,ncol)} directly
  constructs an \dQuote{empty} \linkS4class{sparseMatrix}, as does
  \code{Matrix(FALSE, *)}.

  Although it is sometime possible to mix unclassed matrices (created
  with \code{matrix}) with ones of class \code{"Matrix"}, it is much
  safer to always use carefully constructed ones of class
  \code{"Matrix"}.
}
\seealso{
  The classes \code{\linkS4class{Matrix}},
  \code{\linkS4class{symmetricMatrix}},
  \code{\linkS4class{triangularMatrix}}, and
  \code{\linkS4class{diagonalMatrix}}; further,
  \code{\link{matrix}}.
}
\examples{
Matrix(0, 3, 2)             # 3 by 2 matrix of zeros -> sparse
Matrix(0, 3, 2, sparse=FALSE)# forced 'dense'
Matrix(1:6, 3, 2)           # a 3 by 2 matrix
Matrix(1:6, nrow=3)
Matrix(1:6, ncol=2)

## logical ones:
Matrix(diag(4) >  0)# -> "ldiMatrix" with diag = "U"
Matrix(diag(4) >= 0)# -> "lsyMatrix"
l3 <- upper.tri(matrix(,3,3))
Matrix(l3)  # -> "ltCMatrix"
Matrix(! l3)# -> "ltrMatrix"

Matrix(1:9, nrow=3,
       dimnames = list(c("a", "b", "c"), c("A", "B", "C")))
(I3 <- Matrix(diag(3)))# identity, i.e., unit "diagonalMatrix"
str(I3) # note the empty 'x' slot

(A <- cbind(a=c(2,1), b=1:2))# symmetric *apart* from dimnames
Matrix(A)                    # hence 'dgeMatrix'
(As <- Matrix(A, dimnames = list(NULL,NULL)))# -> symmetric
stopifnot(is(As, "symmetricMatrix"))
}
\keyword{array}
\keyword{algebra}
