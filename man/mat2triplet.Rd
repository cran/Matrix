\name{mat2triplet}
\alias{mat2triplet}
\title{Map Matrix to its Triplet Representation}
\description{
  From an \R object coercible to \code{"\linkS4class{TsparseMatrix}"},
  typically a (sparse) matrix, produce its triplet representation which may
  collapse to a \dQuote{Duplet} in the case of binary aka pattern, such as
  \code{"\linkS4class{nMatrix}"} objects.
}
\usage{
mat2triplet(x, uniqT = FALSE)
}
\arguments{
  \item{x}{any \R object for which \code{as(x, "\linkS4class{TsparseMatrix}")}
    works; typically a \code{\link{matrix}} of one of the \pkg{Matrix}
    package matrices.}
  \item{uniqT}{\code{\link{logical}} indicating if the triplet
    representation should be \sQuote{unique} in the sense of
    \code{\link{uniqTsparse}()}.}
}
\value{
  A \code{\link{list}}, typically with three components,
  \item{i}{vector of   row   indices for all non-zero entries of \code{x}}
  \item{i}{vector of columns indices for all non-zero entries of \code{x}}
  \item{x}{vector of all non-zero entries of \code{x}; exists \bold{only}
    when \code{as(x, "TsparseMatrix")} is \bold{not} a
    \code{"\linkS4class{nsparseMatrix}"}.}

  Note that the \code{\link{order}} of the entries is determined by the
  coercion to \code{"\linkS4class{TsparseMatrix}"} and hence typically
  with increasing \code{j} (and increasing \code{i} within ties of \code{j}).
}
\note{
  The \code{mat2triplet()} utility was created to be a more efficient and
  more predictable substitute for \code{\link{summary}(<sparseMatrix>)}.
  UseRs have wrongly expected the latter to return a data frame with
  columns \code{i} and \code{j} which however is wrong for a
  \code{"\linkS4class{diagonalMatrix}"}.
}
\seealso{
  The \code{summary()} method for \code{"sparseMatrix"},
  \code{\link{summary,sparseMatrix-method}}.

  \code{mat2triplet()} is conceptually the \emph{inverse} function of
  \code{\link{spMatrix}} and (one case of) \code{\link{sparseMatrix}}.
}
\examples{% ../R/sparseMatrix.R
if(FALSE) ## The function is defined (don't redefine here!), simply as
mat2triplet <- function(x, uniqT = FALSE) {
    T <- as(x, "TsparseMatrix")
    if(uniqT && anyDuplicatedT(T)) T <- .uniqTsparse(T)
    if(is(T, "nsparseMatrix"))
         list(i = T@i + 1L, j = T@j + 1L)
    else list(i = T@i + 1L, j = T@j + 1L, x = T@x)
}

i <- c(1,3:8); j <- c(2,9,6:10); x <- 7 * (1:7)
(Ax <- sparseMatrix(i, j, x = x)) ##  8 x 10 "dgCMatrix"
str(trA <- mat2triplet(Ax))
stopifnot(i == sort(trA$i),  sort(j) == trA$j,  x == sort(trA$x))

D <- Diagonal(x=4:2)
summary(D)
str(mat2triplet(D))
}
\keyword{classes}
\keyword{manip}
\keyword{utilities}

