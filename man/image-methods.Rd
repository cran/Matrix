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
\alias{image,lsparseMatrix-method}
\alias{image,nsparseMatrix-method}
\description{
  Methods for function \code{\link[graphics]{image}} in package
  \pkg{Matrix}.  An image of a matrix simply color codes all matrix
  entries and draws the \eqn{n\times m}{n x m} matrix using an
  \eqn{n\times m}{n x m} grid of (colored) rectangles.

  The \pkg{Matrix} package \code{image} methods are based on
  \code{\link[lattice]{levelplot}()} from package \pkg{lattice}; hence
  these methods return an \dQuote{object} of class \code{"trellis"},
  producing a graphic when (auto-) \code{\link{print}()}ed.
}
\usage{% want \usage{} since we have many "surprising arguments"
\S4method{image}{dgTMatrix}(x,
      xlim = c(1, di[2]),
      ylim = c(di[1], 1), aspect = "iso",
      sub = sprintf("Dimensions: \%d x \%d", di[1], di[2]),
      xlab = "Column", ylab = "Row", cuts = 15,
      useRaster = FALSE,
      useAbs = NULL, colorkey = !useAbs,
      col.regions = NULL,
      lwd = NULL, border.col = NULL, \dots)
}
\arguments{
  \item{x}{a Matrix object, i.e., fulfilling \code{\link{is}(x, "Matrix")}.}
  \item{xlim, ylim}{x- and y-axis limits; may be used to \dQuote{zoom
      into} matrix.  Note that \eqn{x,y} \dQuote{feel reversed}:
    \code{ylim} is for the rows (= 1st index) and \code{xlim} for the
    columns (= 2nd index).  For convenience, when the limits are integer
    valued, they are both extended by \code{0.5}; also, \code{ylim} is
    always used decreasingly.}
  \item{aspect}{aspect ratio specified as number (y/x) or string;
    see \code{\link[lattice]{levelplot}}.}
  \item{sub, xlab, ylab}{axis annotation with sensible defaults;
    see \code{\link{plot.default}}.}
  \item{cuts}{number of levels the range of matrix values would be
    divided into.}
  \item{useRaster}{logical indicating if raster graphics should be used
    (instead of the tradition rectangle vector drawing).  If true,
    \code{\link[lattice]{panel.levelplot.raster}} (from \pkg{lattice}
    package) is used, and the colorkey is also done via rasters, see
    also \code{\link[lattice]{levelplot}} and possibly
    \code{\link[grid]{grid.raster}}.

    Note that using raster graphics may often be faster, but can be slower,
    depending on the matrix dimensions and the graphics device (dimensions).}
  \item{useAbs}{logical indicating if \code{\link{abs}(x)} should be
    shown; if \code{TRUE}, the former (implicit) default, the default
    \code{col.regions} will be \code{\link{grey}} colors (and no
    \code{colorkey} drawn).  The default is \code{FALSE} unless the
    matrix has no negative entries.}
  \item{colorkey}{logical indicating if a color key aka \sQuote{legend}
    should be produced.  Default is to draw one, unless \code{useAbs} is
    true.  You can also specify a \code{\link{list}}, see
    \code{\link[lattice]{levelplot}}, such as\code{list(raster=TRUE)} in
    the case of rastering.}
  \item{col.regions}{vector of gradually varying colors; see
    \code{\link[lattice]{levelplot}}.}
  \item{lwd}{(only used when \code{useRaster} is false:) non-negative
    number or \code{NULL} (default), specifying the line-width of the
    rectangles of each non-zero matrix entry (drawn by
    \code{\link[grid]{grid.rect}}).  The default depends on the matrix
    dimension and the device size.}
  \item{border.col}{color for the border of each rectangle. \code{NA}
    means no border is drawn.   When \code{NULL} as by default,
    \code{border.col <- if(lwd < .01) NA else NULL} is used.
    Consider using an opaque color instead of \code{NULL} which
    corresponds to \code{grid::\link[grid]{get.gpar}("col")}.}
  \item{\dots}{further arguments passed to methods and
    \code{\link[lattice]{levelplot}}, notably \code{at} for specifying
    (possibly non equidistant) cut values for dividing the matrix
    values (superseding \code{cuts} above).}% FIXME? example *using*  at=..
}
\section{Methods}{
  All methods currently end up calling the method for the
  \code{\linkS4class{dgTMatrix}} class.
  Use \code{showMethods(image)} to list them all.
}
\value{
  as all \pkg{lattice} graphics functions, \code{image(<Matrix>)}
  returns a \code{"trellis"} object, effectively the result of
  \code{\link[lattice]{levelplot}()}.
}
\seealso{
  \code{\link[lattice]{levelplot}}, and
  \code{\link[lattice]{print.trellis}} from package \pkg{lattice}.
}
\examples{
showMethods(image)
## If you want to see all the methods' implementations:
showMethods(image, incl=TRUE, inherit=FALSE)
\dontshow{## warnings should not happen here, notably when print(<trellis>)
op <- options(warn = 2)
}
data(CAex)
image(CAex, main = "image(CAex)")
image(CAex, useAbs=TRUE, main = "image(CAex, useAbs=TRUE)")

cCA <- Cholesky(crossprod(CAex), Imult = .01)
## See  ?print.trellis --- place two image() plots side by side:
print(image(cCA, main="Cholesky(crossprod(CAex), Imult = .01)"),
      split=c(x=1,y=1,nx=2, ny=1), more=TRUE)
print(image(cCA, useAbs=TRUE),
      split=c(x=2,y=1,nx=2,ny=1))

data(USCounties)
image(USCounties)# huge
image(sign(USCounties))## just the pattern
    # how the result looks, may depend heavily on
    # the device, screen resolution, antialiasing etc
    # e.g. x11(type="Xlib") may show very differently than cairo-based

## Drawing borders around each rectangle;
    # again, viewing depends very much on the device:
image(USCounties[1:400,1:200], lwd=.1)
## Using (xlim,ylim) has advantage : matrix dimension and (col/row) indices:
image(USCounties, c(1,200), c(1,400), lwd=.1)
image(USCounties, c(1,300), c(1,200), lwd=.5 )
image(USCounties, c(1,300), c(1,200), lwd=.01)
## These 3 are all equivalent :
(I1 <- image(USCounties, c(1,100), c(1,100), useAbs=FALSE))
 I2 <- image(USCounties, c(1,100), c(1,100), useAbs=FALSE,        border.col=NA)
 I3 <- image(USCounties, c(1,100), c(1,100), useAbs=FALSE, lwd=2, border.col=NA)
stopifnot(all.equal(I1, I2), all.equal(I2, I3))
## using an opaque border color
image(USCounties, c(1,100), c(1,100), useAbs=FALSE, lwd=3, border.col = adjustcolor("skyblue", 1/2))
\dontshow{options(op)}
if(doExtras <- interactive() || nzchar(Sys.getenv("R_MATRIX_CHECK_EXTRA")) ||
    identical("true", unname(Sys.getenv("R_PKG_CHECKING_doExtras")))) {
## Using raster graphics: For PDF this would give a 77 MB file,
## however, for such a large matrix, this is typically considerably
## *slower* (than vector graphics rectangles) in most cases :
if(doPNG <- !dev.interactive())
   png("image-USCounties-raster.png", width=3200, height=3200)
image(USCounties, useRaster = TRUE) # should not suffer from anti-aliasing
if(doPNG)
   dev.off()
   ## and now look at the *.png image in a viewer you can easily zoom in and out
}#only if(doExtras)
}
\keyword{methods}
\keyword{hplot}
