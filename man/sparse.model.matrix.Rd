\name{sparse.model.matrix}
\title{Construct Sparse Design / Model Matrices}
\alias{sparse.model.matrix}
\alias{fac2sparse}
\alias{fac2Sparse}
\description{Construct a sparse model or \dQuote{design} matrix,
  form a formula and data frame (\code{sparse.model.matrix}) or a single
  factor (\code{fac2sparse}).

  The \code{fac2[Ss]parse()} functions are utilities, also used
  internally in the principal user level function
  \code{sparse.model.matrix()}.
}
\usage{
sparse.model.matrix(object, data = environment(object),
		    contrasts.arg = NULL, xlev = NULL, transpose = FALSE,
		    drop.unused.levels = FALSE, row.names = TRUE,
		    verbose = FALSE, \dots)

fac2sparse(from, to = c("d", "i", "l", "n", "z"),
	   drop.unused.levels = TRUE, giveCsparse = TRUE)
fac2Sparse(from, to = c("d", "i", "l", "n", "z"),
	   drop.unused.levels = TRUE, giveCsparse = TRUE,
	   factorPatt12, contrasts.arg = NULL)
}
\arguments{
  \item{object}{an object of an appropriate class.  For the default
    method, a model formula or terms object.}
  \item{data}{a data frame created with \code{\link{model.frame}}.  If
    another sort of object, \code{model.frame} is called first.}
  \item{contrasts.arg}{\describe{
      \item{for \code{sparse.model.matrix()}:}{A list, whose entries are
	contrasts suitable for input to the \code{\link{contrasts}}
	replacement function and whose names are the names of columns
	of \code{data} containing \code{\link{factor}}s.}
      \item{for \code{fac2Sparse()}:}{character string or \code{NULL} or
	(coercable to) \code{"\linkS4class{sparseMatrix}"}, specifying the
	contrasts to be applied to the factor levels.}
  }}
  \item{xlev}{to be used as argument of \code{\link{model.frame}} if
    \code{data} has no \code{"terms"} attribute.}
  \item{transpose}{logical indicating if the \emph{transpose} should be
    returned; if the transposed is used anyway, setting \code{transpose = TRUE}
    is more efficient.}
  \item{drop.unused.levels}{should factors have unused levels dropped?
    The default for \code{sparse.model.matrix} has been changed to
    \code{FALSE}, 2010-07, for compatibility with \R's standard (dense)
    \code{\link{model.matrix}()}.}
  \item{row.names}{logical indicating if row names should be used.}
  \item{verbose}{logical or integer indicating if (and how much)
    progress output should be printed.}
  \item{\dots}{further arguments passed to or from other methods.}
  \item{from}{(for \code{fac2sparse()}:) a \code{\link{factor}}.}
  \item{to}{a character indicating the \dQuote{kind} of sparse matrix to
    be returned.  The default, \code{"d"} is for \code{\link{double}}.}
  \item{giveCsparse}{(for \code{fac2sparse()}:) logical indicating if
    the result must be a \code{\linkS4class{CsparseMatrix}}.}
  \item{factorPatt12}{logical vector, say \code{fp}, of length two;
    when \code{fp[1]} is true, return \dQuote{contrasted} \code{t(X)};
    when \code{fp[2]} is true, the original (\dQuote{dummy})
    \code{t(X)}, i.e, the result of \code{\link{fac2sparse}()}.}
}
\value{
  a sparse matrix, extending \code{\linkS4class{CsparseMatrix}} (for
  \code{fac2sparse()} if \code{giveCsparse} is true as per default; a
  \code{\linkS4class{TsparseMatrix}}, otherwise).

  For \code{fac2Sparse()}, a \code{\link{list}} of length two, both
  components with the corresponding transposed model matrix, where the
  corresponding \code{factorPatt12} is true.

  Note that \code{\link[MatrixModels]{model.Matrix}(*, sparse=TRUE)}
  from package \pkg{MatrixModels} may be often be preferable to
  \code{sparse.model.matrix()} nowadays, as \code{model.Matrix()}
  returns \code{\link[MatrixModels:modelMatrix-class]{modelMatrix}}
  objects with additional slots \code{assign} and \code{contrasts} which
  relate back to the variables used.

  \code{fac2sparse()}, the basic workhorse of
  \code{sparse.model.matrix()}, returns the \emph{transpose}
  (\code{\link{t}}) of the model matrix.
}
\author{Doug Bates and Martin Maechler, with initial suggestions from Tim
  Hesterberg.
}
\seealso{
  \code{\link{model.matrix}} in standard \R's package \pkg{stats}.\cr
  \code{\link[MatrixModels]{model.Matrix}} which calls
  \code{sparse.model.matrix} or \code{model.matrix} depending on its
  \code{sparse} argument may be preferred to \code{sparse.model.matrix}.

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

(ff <- gl(3,4,, c("X","Y", "Z")))
fac2sparse(ff) #  3 x 12 sparse Matrix of class "dgCMatrix"
##
##  X  1 1 1 1 . . . . . . . .
##  Y  . . . . 1 1 1 1 . . . .
##  Z  . . . . . . . . 1 1 1 1

## can also be computed via sparse.model.matrix():
f30 <- gl(3,0    )
f12 <- gl(3,0, 12)
stopifnot(
  all.equal(t( fac2sparse(ff) ),
	    sparse.model.matrix(~ 0+ff),
	    tolerance = 0, check.attributes=FALSE),
  is(M <- fac2sparse(f30, drop= TRUE),"CsparseMatrix"), dim(M) == c(0, 0),
  is(M <- fac2sparse(f30, drop=FALSE),"CsparseMatrix"), dim(M) == c(3, 0),
  is(M <- fac2sparse(f12, drop= TRUE),"CsparseMatrix"), dim(M) == c(0,12),
  is(M <- fac2sparse(f12, drop=FALSE),"CsparseMatrix"), dim(M) == c(3,12)
 )
}
\keyword{models}
