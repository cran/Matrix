## METHODS FOR GENERIC: image
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## NB: currently going via dgTMatrix in all cases,
##     which is not so inefficient, as levelplot()
##     needs 'i' and 'j' anyway

.image.dgT <-
function(x,
         xlim = c(1, di[2L]),
         ylim = c(di[1L], 1),
         aspect = "iso",
         sub = sprintf("Dimensions: %d x %d", di[1L], di[2L]),
         xlab = "Column",
         ylab = "Row",
         cuts = 15,
         useRaster = FALSE,
         useAbs = NULL,
         colorkey = !useAbs,
         col.regions = NULL,
         lwd = NULL,
         border.col = NULL,
         ...)
{
    ## 'at' can remain missing and be passed to levelplot
    di <- x@Dim
    xx <- x@x
    empty.x <- length(xx) == 0L && length(x) > 0L
    if(empty.x) { # workaround having "empty" matrix
        xx <- 0
        x@i <- x@j <- 0L
    }
    if(missing(useAbs)) ## use abs() when all values are non-neg
        useAbs <- if(empty.x) FALSE else min(xx, na.rm = TRUE) >= 0
    else if(useAbs)
        xx <- abs(xx)
    ## rx <- range(xx, finite = TRUE)
    ## FIXME: make use of 'cuts' now
    ##	    and call levelplot() with 'at = ',
    ##      making sure 0 is included and matching
    ##	    *exactly* - rather than approximately
    if(is.null(col.regions)) {
        l.col <- empty.x || diff(rx <- range(xx, finite = TRUE)) == 0
        col.regions <-
            if(useAbs) {
                grey(if(l.col)
                         0.9
                     else seq(from = 0.7, to = 0, length.out = 100L))
            } else if(l.col)
                "gray90"
            else { ## no abs(.), rx[1] < 0 typically
                nn <- 100
                n0 <- min(nn, max(0, round((0 - rx[1L])/(rx[2L]-rx[1L]) * nn)))
                col.regions <-
                    c(colorRampPalette(c("blue3", "gray80"))(n0),
                      colorRampPalette(c("gray75","red3"))(nn - n0))
            }
    }
    if(!is.null(lwd) && !(is.numeric(lwd) && all(lwd >= 0))) # allow lwd=0
        stop("'lwd' must be NULL or non-negative numeric")
    stopifnot(length(xlim) == 2L, length(ylim) == 2L)
    ## ylim: the rows count from top to bottom:
    ylim <- sort(ylim, decreasing = TRUE)
    if(all(xlim == round(xlim))) xlim <- xlim + c(-0.5, +0.5)
    if(all(ylim == round(ylim))) ylim <- ylim + c(+0.5, -0.5) # decreasing!

    panel <-
    if(useRaster)
        panel.levelplot.raster
    else {
        function(x, y, z, subscripts, at, ..., col.regions) {
            x <- as.numeric(x[subscripts])
            y <- as.numeric(y[subscripts])
            ##
            ## FIXME: use  level.colors() here and 'at' from above --
            ## -----  look at 'zcol' in  panel.levelplot()
            numcol <- length(at) - 1L
            num.r <- length(col.regions)
            col.regions <-
                if(num.r <= numcol)
                    rep_len(col.regions, numcol)
                else col.regions[1 + ((1:numcol-1)*(num.r-1)) %/% (numcol-1)]
            zcol <- rep.int(NA_integer_, length(z))
            for(i in seq_along(col.regions))
                zcol[!is.na(x) & !is.na(y) & !is.na(z) &
                     at[i] <= z & z < at[i+1L]] <- i
            zcol <- zcol[subscripts]

            if(any(subscripts)) {
                ## the line-width used in grid.rect() inside
                ## levelplot()'s panel for the *border* of the
                ## rectangles: levelplot()panel has lwd= 0.01:

                ## Here: use "smart" default !
                if(is.null(lwd)) {
                    wh <- current.viewport()[c("width", "height")]
                    ## wh : current viewport dimension in pixel
                    wh <- (par("cra") / par("cin")) *
                        c(convertWidth(wh$width, "inches",
                                       valueOnly = TRUE),
                          convertHeight(wh$height, "inches",
                                        valueOnly = TRUE))

                    pSize <- wh/di ## size of one matrix-entry in pixels
                    pA <- prod(pSize) # the "area"
                    p1 <- min(pSize)
                    lwd <- ## crude for now
                        if(p1 < 2 || pA < 6) 0.01 # effectively 0
                        else if(p1 >= 4) 1
                        else if(p1 > 3) 0.5
                        else 0.2
                    ## browser()
                    Matrix.msg("rectangle size ",
                               paste(round(pSize, 1L), collapse = " x "),
                               " [pixels];  --> lwd :", formatC(lwd))
                } else stopifnot(is.numeric(lwd), all(lwd >= 0)) # allow 0
                if(is.null(border.col) && lwd < 0.01) # no border
                    border.col <- NA
                grid.rect(x = x, y = y, width = 1, height = 1,
                          default.units = "native",
                          ## FIXME?: allow 'gp' to be passed via '...' !!
                          gp = gpar(fill = col.regions[zcol],
                                    lwd = lwd, col = border.col))
            }
        }
    } # panel <- if(useRaster) ...

    levelplot(xx ~ (x@j + 1L) * (x@i + 1L), # no 'data'
              sub = sub,
              xlab = xlab,
              ylab = ylab,
              xlim = xlim,
              ylim = ylim,
              aspect = aspect,
              colorkey = colorkey,
              col.regions = col.regions,
              cuts = cuts,
              par.settings = list(background = list(col = "transparent")),
              panel = panel,
              ...)
}

setMethod("image", "dgTMatrix", .image.dgT)

setMethod("image", "Matrix",
          function(x, ...)
              image(..sparse2d(.sparse2g(as(x, "TsparseMatrix"))), ...))

setMethod("image", "CHMfactor",
          function(x, ...)
              image(.sparse2g(as(x, "TsparseMatrix")), ...))

rm(.image.dgT)
