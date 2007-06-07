\name{colSums}
\alias{colSums}
\alias{colMeans}
\alias{rowSums}
\alias{rowMeans}
\alias{colMeans,sparseMatrix-method}
\alias{colSums,sparseMatrix-method}
\alias{rowMeans,sparseMatrix-method}
\alias{rowSums,sparseMatrix-method}
%
\title{Form Row and Column Sums and Means}
% see also ~/R/D/r-devel/R/src/library/base/man/colSums.Rd
\description{
  Form row and column sums and means for \code{\linkS4class{Matrix}} objects.
}
\usage{
colSums (x, na.rm = FALSE, dims = 1, \dots)
rowSums (x, na.rm = FALSE, dims = 1, \dots)
colMeans(x, na.rm = FALSE, dims = 1, \dots)
rowMeans(x, na.rm = FALSE, dims = 1, \dots)

\S4method{colSums}{sparseMatrix}(x, na.rm = FALSE,
        dims = 1, sparseResult = FALSE, \dots)
\S4method{rowSums}{sparseMatrix}(x, na.rm = FALSE,
        dims = 1, sparseResult = FALSE, \dots)
\S4method{colMeans}{sparseMatrix}(x, na.rm = FALSE,
        dims = 1, sparseResult = FALSE, \dots)
\S4method{rowMeans}{sparseMatrix}(x, na.rm = FALSE,
        dims = 1, sparseResult = FALSE, \dots)
}
\arguments{
  \item{x}{a Matrix, i.e., inheriting from  \code{\linkS4class{Matrix}}.}
  \item{na.rm}{logical. Should missing values (including \code{NaN})
    be omitted from the calculations?}
  \item{dims}{completely ignored by the \code{Matrix} methods.}
  \item{\dots}{potentially further arguments, for method \code{<->}
    generic compatibility.}
  \item{sparseResult}{logical indicating if the result should be sparse,
    i.e., inheriting from class \code{\linkS4class{sparseVector}}.}
}
% \details{
%   ~~ If necessary, more details than the description above ~~
% }
\value{
  returns a numeric vector if \code{sparseResult} is \code{FALSE} as per
  default.  Otherwise, returns a \code{\linkS4class{sparseVector}}.
}
%\author{Martin}
\seealso{\code{\link[base]{colSums}} and the
  \code{\linkS4class{sparseVector}} classes.
}
\examples{
(M <- bdiag(Diagonal(2), matrix(1:3, 3,4), diag(3:2)))
colSums(M)
d <- Diagonal(10, c(0,0,10,0,2,rep(0,5)))
MM <- kronecker(d, M)
dim(MM) # 40 140
  cm <- colSums(MM)
(scm <- colSums(MM, sparseResult = TRUE))
stopifnot(is(scm, "sparseVector"),
          identical(cm, as.numeric(scm)))
 rowSums(MM, sparseResult = TRUE) # 8 of 40 are not zero
if(FALSE) ## FIXME! -- requires <sparseVec> / <scalar>
rowMeans(MM, sparseResult = TRUE)
%%  rowMeans(MM, sparseResult = TRUE, na.rm = TRUE) % needs even more
}
\keyword{array}
\keyword{algebra}
\keyword{arith}
