\name{is.null.DN}
\title{Are the Dimnames \code{dn} NULL-like ?}
%
\keyword{array}
\keyword{attribute}
\keyword{programming}
\keyword{utilities}
%
\alias{is.null.DN}
%
\description{
  Are the \code{\link{dimnames}} \code{dn} \code{\link{NULL}}-like?

  \code{is.null.DN(dn)} is less strict than \code{\link{is.null}(dn)},
  because it is also true (\code{\link{TRUE}}) when the dimnames
  \code{dn} are \dQuote{like} \code{NULL}, or \code{list(NULL,NULL)}, as
  they can easily be for the traditional \R matrices
  (\code{\link{matrix}}) which have no formal \code{\link{class}}
  definition, and hence much freedom in how their \code{\link{dimnames}}
  look like.
}
\usage{
  is.null.DN(dn)
}
\arguments{
  \item{dn}{\code{\link{dimnames}()} of a \code{\link{matrix}}-like \R
    object.
  }
}
\note{
  This function is really to be used on \dQuote{traditional} matrices
  rather than those inheriting from \code{\linkS4class{Matrix}}, as
  the latter will always have dimnames \code{list(NULL,NULL)} exactly,
  in such a case.
}
\value{
  \code{\link{logical}} \code{\link{TRUE}} or \code{\link{FALSE}}.
}
%% \details{
%% }
\author{Martin Maechler}
\seealso{
  \code{\link{is.null}},
  \code{\link{dimnames}}, \code{\link{matrix}}.
}
\examples{
\dontshow{ % for R_DEFAULT_PACKAGES=NULL
library(stats, pos = "package:base", verbose = FALSE)
library(utils, pos = "package:base", verbose = FALSE)
}
m1 <- m2 <- m3 <- m4 <- m <-
    matrix(round(100 * rnorm(6)), 2, 3)
dimnames(m1) <- list(NULL, NULL)
dimnames(m2) <- list(NULL, character())
dimnames(m3) <- rev(dimnames(m2))
dimnames(m4) <- rep(list(character()),2)

m4 # prints absolutely identically to m

c.o <- capture.output
cm <- c.o(m)
stopifnot(exprs = {
    m == m1; m == m2; m == m3; m == m4
	identical(cm, c.o(m1));	identical(cm, c.o(m2))
	identical(cm, c.o(m3)); identical(cm, c.o(m4))
})

hasNoDimnames <- function(.) is.null.DN(dimnames(.))
stopifnot(exprs = {
    hasNoDimnames(m)
    hasNoDimnames(m1); hasNoDimnames(m2)
    hasNoDimnames(m3); hasNoDimnames(m4)
    hasNoDimnames(Matrix(m) -> M)
    hasNoDimnames(as(M, "sparseMatrix"))
})
}
