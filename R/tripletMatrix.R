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
                   xlim = c(-1, matdim[2]),
                   ylim = c(matdim[1], -1),
                   sub = sprintf("Dimensions: %d x %d", matdim[1], matdim[2]),
                   xlab = "Column", ylab = "Row",
                   cuts = 20,
                   col.regions = grey(seq(from = 0.7, to = 0, length = 100)),
                   ...)
      {
          ##require("lattice", character = TRUE, quietly = TRUE)
          
          matdim <- x@Dim
          lattice::levelplot(abs(x@x) ~ x@j * x@i,
                             sub = sub,
                             xlab = xlab, ylab = ylab,
                             xlim = xlim, ylim = ylim,
                             col.regions = col.regions,
                             par.settings = list(background = list(col = "transparent")),
                             panel = function(x, y, z, zcol, subscripts, ..., col.regions) {
                                 x <- as.numeric(x[subscripts])
                                 y <- as.numeric(y[subscripts])
                                 zcol <- as.numeric(zcol[subscripts])
                                 if (any(subscripts))
                                     grid::grid.rect(x = x, y = y, width = 1, height = 1, 
                                                     default.units = "native",
                                                     gp = grid::gpar(fill = col.regions[zcol], col = NULL))
                             },
                             ...)
      })

























