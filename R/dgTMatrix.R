setAs("dgTMatrix", "dgCMatrix",
      function(from) .Call("dgTMatrix_to_dgCMatrix", from) )

setAs("dgTMatrix", "dgeMatrix",
      function(from) .Call("dgTMatrix_to_dgeMatrix", from) )

setAs("dgTMatrix", "matrix",
      function(from) .Call("dgTMatrix_to_matrix", from) )

setMethod("crossprod", signature(x = "dgTMatrix", y = "missing"),
          function(x, y = NULL)
          .Call("csc_crossprod", as(x, "dgCMatrix")))

setMethod("crossprod", signature(x = "dgTMatrix", y = "matrix"),
          function(x, y = NULL)
          .Call("csc_matrix_crossprod", as(x, "dgCMatrix"), y))

setMethod("crossprod", signature(x = "dgTMatrix", y = "numeric"),
          function(x, y = NULL)
          .Call("csc_matrix_crossprod", as(x, "dgCMatrix"), as.matrix(y)))

setMethod("tcrossprod", signature(x = "dgTMatrix"),
          function(x)
          .Call("csc_tcrossprod", as(x, "dgCMatrix")))

setMethod("image", "dgTMatrix",
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

## Uses the triplet convention of *adding* entries with same (i,j):
setMethod("+", signature(e1 = "dgTMatrix", e2 = "dgTMatrix"),
          function(e1, e2) {
              if (any(e1@Dim != e2@Dim))
                  stop("Dimensions not compatible for addition")
              new("dgTMatrix", i = c(e1@i, e2@i), j = c(e1@j, e2@j),
                  x = c(e1@x, e2@x), Dim = e1@Dim)
          })

setMethod("t", signature(x = "dgTMatrix"),
          function(x)
          new("dgTMatrix", i = x@j, j = x@i, x = x@x, Dim = rev(x@Dim)))

setMethod("isSymmetric", signature(object = "dgTMatrix"),
          ## This is not a complete test.  Probably use .Call for complete test.
          function(object, ...) {
              i <- object@i
              j <- object@j
              all(sort(paste(i, j, sep=':')) == sort(paste(j, i, sep=':')))
          })

setAs("dgTMatrix", "dsCMatrix",
      function(from) {
          i <- from@i
          j <- from@j
          if (any(upper <- j > i)) {
              from@j[upper] <- i[upper]
              from@i[upper] <- j[upper]
          }
          as(as(from, "dgCMatrix"), "dsCMatrix")
      })

