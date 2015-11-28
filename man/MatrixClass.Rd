\name{MatrixClass}
\title{The Matrix (Super-) Class of a Class}
\alias{MatrixClass}
\description{
  Return the (maybe super-)\code{\link{class}} of class \code{cl} from
  package \pkg{Matrix}, returning \code{\link{character}(0)} if there is none.
}
\usage{
MatrixClass(cl, cld = getClassDef(cl), ...Matrix = TRUE,
            dropVirtual = TRUE, ...)
}
\arguments{
  \item{cl}{string, class name}
  \item{cld}{its class definition}
  \item{...Matrix}{\code{\link{logical}} indicating if the result must be of
    pattern \code{"[dlniz]..Matrix"} where the first letter "[dlniz]"
    denotes the content kind.}
  \item{dropVirtual}{\code{\link{logical}} indicating if virtual classes
    are included or not.}% ?? (FIXME) -- example
  \item{\dots}{further arguments are passed to
    \code{\link{.selectSuperClasses}()}.}
}
\value{
  a \code{\link{character}} string
}
\author{Martin Maechler, 24 Mar 2009}
%% \details{
%% }
\seealso{
  \code{\linkS4class{Matrix}}, the mother of all \pkg{Matrix} classes.
}
\examples{
mkA <- setClass("A", contains="dgCMatrix")
(A <- mkA())
stopifnot(identical(
     MatrixClass("A"),
     "dgCMatrix"))
}
\keyword{classes}
