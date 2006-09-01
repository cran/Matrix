## Now in ./Tsparse.R
## setAs("dgTMatrix", "dgCMatrix",
##       function(from) .Call(Tsparse_to_Csparse, from, FALSE)
##       )

setAs("dgTMatrix", "dgeMatrix",
      function(from) .Call(dgTMatrix_to_dgeMatrix, from)
      )

setAs("dgTMatrix", "matrix",
      function(from) .Call(dgTMatrix_to_matrix, from)
      )

setAs("dgeMatrix", "dgTMatrix",
      function(from) as(as(from, "dgCMatrix"), "dgTMatrix"))

setAs("dgTMatrix", "dsCMatrix",
      function(from) {
          if (!isSymmetric(from))
              stop("cannot coerce non-symmetric dgTMatrix to dsCMatrix class")
          upper <- from@i <= from@j
          uC <- as(new("dgTMatrix", Dim = from@Dim,  Dimnames = from@Dimnames,
                       i = from@i[upper],
                       j = from@j[upper], x = from@x[upper]), "dgCMatrix")
          new("dsCMatrix", Dim = uC@Dim, p = uC@p, i = uC@i, x = uC@x, uplo = "U")
      })

setAs("dgTMatrix", "dsTMatrix",
      function(from) {
	  if(isSymmetric(from)) {
	      upper <- from@i <= from@j
	      new("dsTMatrix", Dim = from@Dim, Dimnames = from@Dimnames,
		  i = from@i[upper],
		  j = from@j[upper], x = from@x[upper], uplo = "U")
	  }
	  else stop("not a symmetric matrix")})

setAs("dgTMatrix", "dtTMatrix",
      function(from) {
	  if(isTr <- isTriangular(from))
	      gt2tT(from, uplo = attr(isTr, "kind"),
		    diag = "N") ## improve: also test for unit diagonal
	  else stop("not a triangular matrix")})

setAs("matrix", "dgTMatrix",
      function(from) {
	  x <- as.double(from)
	  nz <- as.logical(x)
	  new("dgTMatrix", Dim = dim(from),
	      i = row(from)[nz] - 1:1,
	      j = col(from)[nz] - 1:1,
	      x = x[nz])
      })



## "[" methods are now in ./Tsparse.R

## "[<-" methods { setReplaceMethod()s }  too ...


## "crossprod" methods too ...
## setMethod("crossprod", signature(x = "dgTMatrix", y = "missing"),
##           function(x, y = NULL)
##           .Call(csc_crossprod, as(x, "dgCMatrix")))

## setMethod("crossprod", signature(x = "dgTMatrix", y = "matrix"),
##           function(x, y = NULL)
##           .Call(csc_matrix_crossprod, as(x, "dgCMatrix"), y))

##setMethod("crossprod", signature(x = "dgTMatrix", y = "numeric"),
##          function(x, y = NULL)
##          .Call(csc_matrix_crossprod, as(x, "dgCMatrix"), as.matrix(y)))

## setMethod("tcrossprod", signature(x = "dgTMatrix", y = "missing"),
##           function(x, y = NULL)
##           .Call(csc_tcrossprod, as(x, "dgCMatrix")))

setMethod("image", "dgTMatrix",
          function(x,
		   xlim = .5 + c(0, matdim[2]),
		   ylim = .5 + c(matdim[1], 0),
                   sub = sprintf("Dimensions: %d x %d", matdim[1], matdim[2]),
                   xlab = "Column", ylab = "Row",
                   cuts = 20,
                   col.regions = grey(seq(from = 0.7, to = 0, length = 100)),
                   ...)
      {
          matdim <- x@Dim
          levelplot(abs(x@x) ~ (x@j + 1:1) * (x@i + 1:1),
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
                    num.r <- length(col.regions)
		    col.regions <-
			if (num.r <= numcol)
			    rep(col.regions, length = numcol)
			else col.regions[1+ ((1:numcol-1)*(num.r-1)) %/% (numcol-1)]
                    zcol <- rep.int(NA, length(z)) #numeric(length(z))
                    for (i in seq(along = col.regions))
                        zcol[!is.na(x) & !is.na(y) & !is.na(z) &
                             at[i] <= z & z < at[i+1]] <- i

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
              dimCheck(e1, e2)
              new("dgTMatrix", i = c(e1@i, e2@i), j = c(e1@j, e2@j),
                  x = c(e1@x, e2@x), Dim = e1@Dim)
          })

setMethod("t", signature(x = "dgTMatrix"),
          function(x)
          new("dgTMatrix", i = x@j, j = x@i, x = x@x, Dim = rev(x@Dim)))

setMethod("kronecker", signature(X = "dgTMatrix", Y = "dgTMatrix"),
          function (X, Y, FUN = "*", make.dimnames = FALSE, ...)
      {
          if (FUN != "*") stop("kronecker method must use default 'FUN'")
          ydim <- Y@Dim
          xi <- X@i
          xnnz <- length(xi)
          yi <- Y@i
          ynnz <- length(yi)
          new("dgTMatrix", Dim = X@Dim * ydim,
              i = rep.int(yi, xnnz) + ydim[1] * rep.int(xi, rep.int(ynnz, xnnz)),
              j = rep.int(Y@j, xnnz) + ydim[2] * rep.int(X@j, rep.int(ynnz, xnnz)),
              x = as.vector(outer(Y@x, X@x)))
      }, valueClass = "dgTMatrix")

setMethod("writeHB", signature(obj = "dgTMatrix"),
          function(obj, file, ...)
          .Call(Matrix_writeHarwellBoeing, obj, as.character(file), "DGT"))

setMethod("writeMM", signature(obj = "dgTMatrix"),
          function(obj, file, ...)
          .Call(Matrix_writeMatrixMarket, obj, as.character(file), "DGT"))


setMethod("colSums", signature(x = "dgTMatrix"),
	  function(x, na.rm = FALSE, dims = 1)
          tapply1(x@x, factor(x@j, 0:(x@Dim[2]-1)), sum, na.rm = na.rm),
	  valueClass = "numeric")
setMethod("colMeans", signature(x = "dgTMatrix"),
	  function(x, na.rm = FALSE, dims = 1)
          tapply1(x@x, factor(x@j, 0:(x@Dim[2]-1)), mean, na.rm = na.rm),
	  valueClass = "numeric")

setMethod("rowSums", signature(x = "dgTMatrix"),
	  function(x, na.rm = FALSE, dims = 1)
          tapply1(x@x, factor(x@i, 0:(x@Dim[1]-1)), sum, na.rm = na.rm),
	  valueClass = "numeric")
setMethod("rowMeans", signature(x = "dgTMatrix"),
	  function(x, na.rm = FALSE, dims = 1)
          tapply1(x@x, factor(x@i, 0:(x@Dim[1]-1)), mean, na.rm = na.rm),
	  valueClass = "numeric")
