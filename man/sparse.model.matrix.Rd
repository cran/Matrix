\name{sparse.model.matrix}
\Rdversion{1.1}
\alias{sparse.model.matrix}
\title{Construct Sparse Design / Model Matrices}
\description{
  Construct a Model or \dQuote{Design} Matrix
}
\usage{
sparse.model.matrix(object, data = environment(object),
                    contrasts.arg = NULL, xlev = NULL, transpose = FALSE, ...)
}
\arguments{
  \item{object}{an object of an appropriate class.  For the default
    method, a model formula or terms object.}
  \item{data}{a data frame created with \code{\link{model.frame}}.  If
    another sort of object, \code{model.frame} is called first.}
  \item{contrasts.arg}{A list, whose entries are contrasts suitable for
    input to the \code{\link{contrasts}} replacement function and
    whose names are the names of columns of \code{data} containing
    \code{\link{factor}}s.}
  \item{xlev}{to be used as argument of \code{\link{model.frame}} if
    \code{data} has no \code{"terms"} attribute.}
  \item{transpose}{logical indicating if the \emph{transpose} should be
    returned; if the transposed is used anyway, setting \code{transpose = TRUE}
    is more efficient.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\value{
  a (special) sparse matrix, namely inheriting from class
  \code{\linkS4class{sparseModelMatrix}} which extends
  \code{\linkS4class{CsparseMatrix}}.

  Note that such \code{\linkS4class{modelMatrix}} objects have
  additional slots \code{assign} and \code{contrasts} which relate back
  to the variables used.
}
\author{Doug Bates and Martin Maechler, with initial suggestions from Tim
  Hesterberg.
}
\seealso{
  \code{\link{model.matrix}} in standard \R's package \pkg{stats}.\cr
  \code{\link{model.Matrix}} which calls \code{sparse.model.matrix} or
  \code{model.matrix} depending on its \code{sparse} argument.

  \code{as(f, "sparseMatrix")} (see \code{coerce(from = "factor", ..)}
  in the class doc \linkS4class{sparseMatrix}) produces the
  \emph{transposed} sparse model matrix for a single factor \code{f}
  (and \emph{no} contrasts).
}
\examples{
dd <- data.frame(a = gl(3,4), b = gl(4,1,12))# balanced 2-way
options("contrasts") # the default:  "contr.treatment"
sparse.model.matrix(~ a + b, dd)
sparse.model.matrix(~ -1+ a + b, dd)# no intercept --> even sparser
sparse.model.matrix(~ a + b, dd, contrasts = list(a="contr.sum"))
sparse.model.matrix(~ a + b, dd, contrasts = list(b="contr.SAS"))

## Sparse method is equivalent to the traditional one :
stopifnot(all(sparse.model.matrix(~ a + b, dd) ==
	      Matrix(model.matrix(~ a + b, dd), sparse=TRUE)),
          all(sparse.model.matrix(~ 0+ a + b, dd) ==
	      Matrix(model.matrix(~ 0+ a + b, dd), sparse=TRUE)))
%% many more and tougher examples ---> ../tests/spModel.matrix.R
}
\keyword{models}

