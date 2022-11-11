\name{dimScale}
\title{Scale the Rows and Columns of a Matrix}
\alias{dimScale}
\alias{rowScale}
\alias{colScale}
\description{
  \code{dimScale}, \code{rowScale}, and \code{colScale} implement
  \code{D1 \%*\% x \%*\% D2}, \code{D \%*\% x}, and \code{x \%*\% D}
  for diagonal matrices \code{D1}, \code{D2}, and \code{D} with
  diagonal entries \code{d1}, \code{d2}, and \code{d}, respectively.
  Unlike the explicit products, these functions preserve \code{dimnames(x)}
  and symmetry where appropriate.
}
\usage{
dimScale(x, d1 = sqrt(1/diag(x, names = FALSE)), d2 = d1)
rowScale(x, d)
colScale(x, d)
}
\arguments{
  \item{x}{a matrix, possibly inheriting from virtual class
    \code{\linkS4class{Matrix}}.}
  \item{d1,d2,d}{numeric vectors giving factors by which to scale
    the rows or columns of \code{x}; they are recycled as necessary.}
}
\details{
  \code{dimScale(x)} (with \code{d1} and \code{d2} unset) is only
  roughly equivalent to \code{\link{cov2cor}(x)}.  \code{cov2cor}
  sets the diagonal entries of the result to 1 (exactly);
  \code{dimScale} does not.
}
\value{
  The result of scaling \code{x}, currently always inheriting from
  virtual class \code{\linkS4class{dMatrix}}.

  It inherits from \code{\linkS4class{triangularMatrix}} if and only
  if \code{x} does.  In the special case of \code{dimScale(x, d1, d2)}
  with identical \code{d1} and \code{d2}, it inherits from
  \code{\linkS4class{symmetricMatrix}} if and only if \code{x} does.
}
\author{Mikael Jagan}
\seealso{\code{\link{cov2cor}}}
\examples{
n <- 6L
(x <- forceSymmetric(matrix(1, n, n)))
dimnames(x) <- rep.int(list(letters[seq_len(n)]), 2L)

d <- seq_len(n)
(D <- Diagonal(x = d))

(scx <- dimScale(x, d)) # symmetry and 'dimnames' kept
(mmx <- D \%*\% x \%*\% D) # symmetry and 'dimnames' lost
stopifnot(identical(unname(as(scx, "generalMatrix")), mmx))

rowScale(x, d)
colScale(x, d)
}
