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
                   xlim = c(0, matdim[2] + 1),
                   ylim = c(matdim[1] + 1, 0),
                   sub = sprintf("Dimensions: %d x %d", matdim[1], matdim[2]),
                   xlab = "Column", ylab = "Row",
                   cuts = 20,
                   col.regions = grey(seq(from = 0.7, to = 0, length = 100)),
                   ...) {
              require("lattice", character = TRUE, quietly = TRUE)
              
              matdim <- x@Dim
              levelplot(abs(x@x) ~ x@j * x@i,
                        sub = sub,
                        xlab = xlab, ylab = ylab,
                        xlim = xlim, ylim = ylim,
                        col.regions = col.regions,
                        par.settings = list(background = list(col = "transparent")),
                        ...)
          })
