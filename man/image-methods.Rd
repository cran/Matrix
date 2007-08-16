\name{image-methods}
\title{Methods for image() in Package 'Matrix'}
\docType{methods}
%\alias{image}% <- needed for \usage{.}
\alias{image-methods}
\alias{image,ANY-method}
\alias{image,CHMfactor-method}
\alias{image,Matrix-method}
\alias{image,dgRMatrix-method}
\alias{image,dgCMatrix-method}
\alias{image,dgTMatrix-method}
\alias{image,dsparseMatrix-method}
\alias{image,lgCMatrix-method}
\alias{image,lsCMatrix-method}
\alias{image,lsTMatrix-method}
\alias{image,ltCMatrix-method}
\alias{image,ltTMatrix-method}
\alias{image,ngCMatrix-method}
\alias{image,nsCMatrix-method}
\alias{image,nsTMatrix-method}
\alias{image,ntCMatrix-method}
\alias{image,ntTMatrix-method}
\description{
  Methods for function \code{\link[graphics]{image}} in package
  \pkg{Matrix}.  An image of a matrix

}
\usage{% want \usage{} since we have many "surprising arguments"
\S4method{image}{dgTMatrix}(x,
      xlim = .5 + c(0, di[2]),
      ylim = .5 + c(di[1], 0), aspect = "iso",
      sub = sprintf("Dimensions: \%d x \%d", di[1], di[2]),
      xlab = "Column", ylab = "Row",
      col.regions = grey(seq(from = 0.7, to = 0, length = 100)),
      colorkey = FALSE, \dots)
}
\arguments{
  \item{x}{a Matrix object, i.e., fulfilling
    \code{\link{is}(x, "Matrix")}.}
  \item{xlim, ylim}{x- and y-axis limits; may be used to \dQuote{zoom
      into} matrix.}
  \item{aspect}{aspect ratio specified as number (y/x) or string;
    see \code{\link[lattice]{levelplot}}.}
  \item{sub, xlab, ylab}{axis annotation with sensible defaults;
    see \code{\link{plot.default}}.}
  \item{col.regions}{vector of gradually varying colors; see
    \code{\link[lattice]{levelplot}}.}
  \item{colorkey}{logical indicating if a color key aka \sQuote{legend}
    should be produced.  Defaults to \code{FALSE}.}
  \item{\dots}{further arguments passed to methods and
    \code{\link[lattice]{levelplot}}.}
}
\section{Methods}{
  All methods currently end up calling the method for the
  \code{\linkS4class{dgTMatrix}} class.
  Use \code{showMethods(image)} to list them all.
}
\examples{
showMethods(image)
## If you want to see all the methods' implementations:
showMethods(image, incl=TRUE, inherit=FALSE)
}
\keyword{methods}
\keyword{hplot}
