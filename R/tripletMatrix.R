setAs("tripletMatrix", "cscMatrix",
      function(from)
      .Call("triplet_to_csc", from, PACKAGE = "Matrix")
      )

setAs("tripletMatrix", "geMatrix",
      function(from)
      .Call("triplet_to_geMatrix", from, PACKAGE = "Matrix")
      )

setMethod("image", "tripletMatrix",
          function(x,
                   xlim = c(-0.5, matdim[2]-0.5),
                   ylim = c(matdim[1]-0.5, -0.5),
                   sub = sprintf("Dimensions: %d x %d", matdim[1], matdim[2]),
                   xlab = "Column", ylab = "Row",
                   cuts = 20,
                   col.regions = grey(seq(from = 0.7, to = 0, length = 100)),
                   ...)
      {
          matdim <- x@Dim
          levelplot(abs(x@x) ~ x@j * x@i,
                    sub = sub,
                    xlab = xlab, ylab = ylab,
                    xlim = xlim, ylim = ylim,
                    col.regions = col.regions,
                    par.settings = list(background = list(col = "transparent")),
                    panel = function(x, y, z, subscripts, at, ..., col.regions)
                {
                    x <- as.numeric(x[subscripts])
                    y <- as.numeric(y[subscripts])
                    
                    numcol <- length(at) - 1
                    numcol.r <- length(col.regions)
                    col.regions <-
                        if (numcol.r <= numcol)
                            rep(col.regions, length = numcol)
                        else col.regions[floor(1+(1:numcol-1)*(numcol.r-1)/
                                               (numcol-1))]
                    zcol <- rep(NA, length(z)) #numeric(length(z))
                    for (i in seq(along = col.regions))
                        zcol[!is.na(x) & !is.na(y) & !is.na(z) &
                             z>=at[i] & z<at[i+1]] <- i
                    
                    zcol <- as.numeric(zcol[subscripts])
                    if (any(subscripts))
                        grid.rect(x = x, y = y, width = 1, height = 1, 
                                  default.units = "native",
                                  gp = gpar(fill = col.regions[zcol],
                                  col = NULL))
                }, ...)
      })



























