\name{cBind}
\alias{cBind}
\alias{rBind}
%no \alias{cbind}
%no \alias{rbind}
%
\title{Versions of 'cbind' and 'rbind' recursively built on cbind2/rbind2}
\description{
  The base functions \code{\link{cbind}} and \code{\link{rbind}} are
  defined for an arbitrary number of arguments and hence have the first
  formal argument \code{...}.  For that reason, in the past S4 methods
  could easily be defined for binding together matrices inheriting from
  \code{\link{Matrix}}.

  For that reason, \code{\link{cbind2}} and \code{\link{rbind2}} have
  been provided for binding together \emph{two} matrices, and we have
  defined methods for these and the \code{'Matrix'}-matrices.

  Before \R version 3.2.0 (April 2015), we have needed a substitute for
  \emph{S4-enabled} versions of \code{cbind} and \code{rbind}, and
  provided \code{cBind} and \code{rBind} with identical syntax and
  semantic in order to bind together multiple matrices (\code{"matrix"}
  or \code{"Matrix"} and vectors.
  With \R version 3.2.0 and newer, \code{cBind} and \code{rBind} are
  \emph{deprecated} and produce a deprecation warning (via
  \code{\link{.Deprecated}}), and your code should start using
  \code{cbind()} and \code{rbind()} instead.

}
\usage{
%no cbind(..., deparse.level = 1)
%no rbind(..., deparse.level = 1)
cBind(..., deparse.level = 1)
rBind(..., deparse.level = 1)
}
\arguments{
  \item{\dots}{matrix-like \R objects to be bound together, see
    \code{\link{cbind}} and \code{\link{rbind}}.}
  \item{deparse.level}{integer determining under which circumstances
    column and row names are built from the actual arguments'
    \sQuote{expression}, see \code{\link{cbind}}.}
}
\details{
  The implementation of these is \emph{recursive}, calling
  \code{\link{cbind2}} or
  \code{\link{rbind2}} respectively, where these have methods
  defined and so should dispatch appropriately.
}
\value{
  typically a \sQuote{matrix-like} object of a similar
  \code{\link{class}} as the first argument in \code{\dots}.

  Note that sometimes by default, the result is a
  \code{\linkS4class{sparseMatrix}} if one of the arguments is (even in
  the case where this is not efficient).  In other cases,
  the result is chosen to be sparse when there are more zero entries is
  than non-zero ones (as the default \code{sparse} in
  \code{\link{Matrix}()}).
}
\author{Martin Maechler}
\seealso{\code{\link{cbind2}}, \code{\link{cbind}},
  \code{\link[methods]{Methods}}.
}
\examples{
(a <- matrix(c(2:1,1:2), 2,2))
D <- Diagonal(2)
if(getRversion() < "3.2.0") {
 M1 <- cbind(0, rBind(a, 7))
 print(M1) # remains traditional matrix

 M2 <- cBind(4, a, D, -1, D, 0) # a sparse Matrix
 print(M2)

} else { ## newer versions of R do not need cBind / rBind:

 M1 <- cbind(0, suppressWarnings(rBind(a, 7)))
 print(M1) # remains traditional matrix

 M2 <- suppressWarnings(cBind(4, a, D, -1, D, 0)) # a sparse Matrix
 print(M2)

 stopifnot(identical(M1, cbind(0, rbind(a, 7))),
           identical(M2, cbind(4, a, D, -1, D, 0)))
}# R >= 3.2.0
}
\keyword{array}
\keyword{manip}

