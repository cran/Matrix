\name{drop0}
\alias{drop0}
\title{Drop "Explicit Zeroes" from a Sparse Matrix}
\description{
  Returns a sparse matrix with no \dQuote{explicit zeroes}, i.e., all
  zero or \code{FALSE} entries are dropped from the explicitly indexed
  matrix entries.
}
\usage{
drop0(x, clx = c(class(x)))
}
\arguments{
  \item{x}{a Matrix, typically sparse, i.e., inheriting from
    \code{\linkS4class{sparseMatrix}}.}
  \item{clx}{[optional:] the \code{\link{class}} or \dQuote{class
      definition} (see \code{\link{getClassDef}}; it is of class
    \code{\linkS4class{classRepresentation}}) of \code{x}.\cr
    This argument just exists for the possibility of speedup.}
}
% \details{
%   ~~ If necessary, more details than the description above ~~
% }
\value{
  a Matrix like \code{x} but with no explicit zeros, i.e.,
  \code{!any(x@x == 0)}, always inheriting from
  \code{\linkS4class{CsparseMatrix}}.
}
% \author{Martin}
\seealso{\code{\link{spMatrix}}, class \code{\linkS4class{sparseMatrix}}
}
\examples{
m <- spMatrix(10,20, i= 1:8, j=2:9,
              x = c(0:2,3:-1))
m
drop0(m)
}
\keyword{utilities}
\keyword{array}

