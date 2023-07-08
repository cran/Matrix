\name{invertPerm}
\title{Utilities for Permutation Vectors}
%
\keyword{utilities}
%
\alias{invertPerm}
\alias{signPerm}
\alias{isPerm}
\alias{asPerm}
\alias{invPerm}
%
\description{
  \code{invertPerm} and \code{signPerm} compute the inverse and sign
  of a length-\code{n} permutation vector.  \code{isPerm} tests
  if a length-\code{n} integer vector is a valid permutation vector.
  \code{asPerm} coerces a length-\code{m} transposition vector to a
  length-\code{n} permutation vector, where \code{m <= n}.
}
\usage{
invertPerm(p, off = 1L, ioff = 1L)
signPerm(p, off = 1L)
isPerm(p, off = 1L)
asPerm(pivot, off = 1L, ioff = 1L, n = length(pivot))

invPerm(p, zero.p = FALSE, zero.res = FALSE)
}
\arguments{
  \item{p}{an integer vector of length \code{n}.}
  \item{pivot}{an integer vector of length \code{m}.}
  \item{off}{an integer offset, indicating that \code{p} is
    a permutation of \code{off+0:(n-1)} or that \code{pivot}
    contains \code{m} values sampled with replacement from
    \code{off+0:(n-1)}.}
  \item{ioff}{an integer offset, indicating that the result
    should be a permutation of \code{ioff+0:(n-1)}.}
  \item{n}{a integer greater than or equal to \code{m},
    indicating the length of the result.  Transpositions
    are applied to a permutation vector vector initialized
    as \code{seq_len(n)}.}
  \item{zero.p}{a logical.  Equivalent to \code{off=0} if \code{TRUE}
    and \code{off=1} if \code{FALSE}.}
  \item{zero.res}{a logical.  Equivalent to \code{ioff=0} if \code{TRUE}
    and \code{ioff=1} if \code{FALSE}.}
}
\details{
  \code{invertPerm(p, off, ioff=1)} is equivalent to
  \code{\link{order}(p)} or \code{\link{sort.list}(p)}
  for all values of \code{off}.  For the default value
  \code{off=1}, it returns the value of \code{p} after
  \code{p[p] <- seq_along(p)}.
  
  \code{invPerm} is a simple wrapper around \code{invertPerm},
  retained for backwards compatibility.
}
\value{
  By default, i.e., with \code{off=1} and \code{ioff=1}:
  
  \code{invertPerm(p)} returns an integer vector of length
  \code{length(p)} such that \code{p[invertPerm(p)]}
  and \code{invertPerm(p)[p]} are both \code{seq_along(p)},
  i.e., the identity permutation.

  \code{signPerm(p)} returns 1 if \code{p} is an even permutation
  and \code{-1} otherwise (i.e., if \code{p} is odd).

  \code{isPerm(p)} returns \code{TRUE} if \code{p} is a
  permutation of \code{seq_along(p)} and \code{FALSE} otherwise.

  \code{asPerm(pivot)} returns the result of transposing elements
  \code{i} and \code{pivot[i]} of a permutation vector initialized
  as \code{seq_len(n)}, for \code{i} in \code{seq_along(pivot)}.
}
\seealso{
  Class \code{\linkS4class{pMatrix}} of permutation matrices.
}
\examples{
p <- sample(10L) # a random permutation vector
ip <- invertPerm(p)
s <- signPerm(p)

## 'p' and 'ip' are indeed inverses:
stopifnot(exprs = {
    isPerm(p)
    isPerm(ip)
    identical(s, 1L) || identical(s, -1L)
    identical(s, signPerm(ip))
    identical(p[ip], 1:10)
    identical(ip[p], 1:10)
    identical(invertPerm(ip), p)
})

## Product of transpositions (1 2)(2 1)(4 3)(6 8)(10 1) = (3 4)(6 8)(1 10)
pivot <- c(2L, 1L, 3L, 3L, 5L, 8L, 7L, 8L, 9L, 1L)
q <- asPerm(pivot)
stopifnot(exprs = {
    identical(q, c(10L, 2L, 4L, 3L, 5L, 8L, 7L, 6L, 9L, 1L))
    identical(q[q], seq_len(10L)) # because the permutation is odd:
    signPerm(q) == -1L
})

invPerm # a less general version of 'invertPerm'
\dontshow{
stopifnot(exprs = {
    identical(isPerm(0L), FALSE)
    identical(signPerm(1:2),  1L)
    identical(signPerm(2:1), -1L)
    identical(invertPerm(c(3, 1:2)), c(2:3, 1L)) # 'p' of type "double",
    tryCatch(invPerm(NA), error = function(e) TRUE) # was a segfault
})
}
}
