\name{sparseVector}
\alias{sparseVector}
\title{Sparse Vector Construction from Nonzero Entries}
\description{
  User friendly construction of sparse vectors,
  i.e., objects inheriting from \code{\link{class}}
  \code{\linkS4class{sparseVector}}, from indices and values of its
  non-zero entries.
}
\details{
  zero entries in \code{x} are dropped automatically, analogously as
  \code{\link{drop0}()} acts on sparse matrices.
}
\usage{
sparseVector(x, i, length)
}
\arguments{
  \item{x}{vector of the non zero entries; may be missing in which case a
    \code{"nsparseVector"} will be returned.}
  \item{i}{integer vector (of the same length as \code{x}) specifying
    the indices of the non-zero (or non-\code{TRUE}) entries of the
    sparse vector.}
  \item{length}{length of the sparse vector.}
}
\value{
  a sparse vector, i.e., inheriting from \code{\link{class}}
  \code{\linkS4class{sparseVector}}.
}
\author{Martin Maechler}
\seealso{
  \code{\link{sparseMatrix}()} constructor for sparse matrices;
  the class \code{\linkS4class{sparseVector}}.
}
\examples{
str(sv <- sparseVector(x = 1:10, i = sample(999, 10), length=1000))

sx <- c(0,0,3, 3.2, 0,0,0,-3:1,0,0,2,0,0,5,0,0)
ss <- as(sx, "sparseVector")
stopifnot(identical(ss,
   sparseVector(x = c(2, -1, -2, 3, 1, -3, 5, 3.2),
                i = c(15L, 10:9, 3L,12L,8L,18L, 4L), length = 20L)))

(ns <- sparseVector(i= c(7, 3, 2), length = 10))
stopifnot(identical(ns,
      new("nsparseVector", length = 10, i = c(2, 3, 7))))
}
\keyword{array}
