setAs("dgTMatrix", "dgCMatrix",
      function(from) .Call("dgTMatrix_to_dgCMatrix", from) )

setAs("dgTMatrix", "dgeMatrix",
      function(from) .Call("dgTMatrix_to_dgeMatrix", from) )

setAs("dgTMatrix", "matrix",
      function(from) .Call("dgTMatrix_to_matrix", from) )

setAs("dgeMatrix", "dgTMatrix",
      function(from) as(as(from, "dgCMatrix"), "dgTMatrix"))

setMethod("[", signature(x = "dgTMatrix",
			 i = "missing", j = "missing", drop = "ANY"),
	  function (x, i, j, ..., drop) x)

setMethod("[", signature(x = "dgTMatrix", i = "numeric", j = "missing",
			 drop = "logical"),
	  function (x, i, j, ..., drop) { ## select rows
	      storage.mode(i) <- "integer"
              xi <- x@i + 1:1 # 1-indexing
	      sel <- xi %in% i
              nd <- c(length(i), ncol(x))
	      x <- new("dgTMatrix", Dim = nd,
		       i = match(xi[sel], i) - 1:1,
		       j = x@j[sel],
		       x = x@x[sel])
	      if (drop && any(nd == 1)) drop(as(x,"matrix")) else x
	  })


setMethod("[", signature(x = "dgTMatrix", i = "missing", j = "numeric",
			 drop = "logical"),
	  function (x, i, j, ..., drop) { ## select columns
	      storage.mode(j) <- "integer"
              xj <- x@j + 1:1 # 1-indexing
	      sel <- xj %in% j
              nd <- c(nrow(x), length(j))
	      x <- new("dgTMatrix", Dim = nd,
		       i = x@i[sel],
		       j = match(xj[sel], j) - 1:1,
		       x = x@x[sel])
	      if (drop && any(nd == 1)) drop(as(x,"matrix")) else x
	  })


## How can we get at   A[ ij ]	where ij is (i,j) 2-column matrix?
##  and                A[ LL ]	where LL is a logical *vector*

setMethod("[", signature(x = "dgTMatrix", i = "numeric", j = "missing",
			 drop = "missing"),
	  function(x,i,j, ..., drop)
	  callGeneric(x, i=i, drop = TRUE)# or rather 'FALSE' ?
	  )
setMethod("[", signature(x = "dgTMatrix", i = "missing", j = "numeric",
			 drop = "missing"),
	  function(x,i,j, ..., drop)
	  callGeneric(x, j=j, drop = TRUE)# or rather 'FALSE' ?
	  )

## [.data.frame has : drop = if (missing(i)) TRUE else length(cols) == 1)

setMethod("[", signature(x = "dgTMatrix",
			 i = "numeric", j = "numeric", drop = "logical"),
	  function (x, i, j, ..., drop)
      {
	  ## (i,j, drop) all specified
	  storage.mode(i) <- "integer"
	  storage.mode(j) <- "integer"
          xi <- x@i + 1:1
          xj <- x@j + 1:1
	  sel <- (xi %in% i) & (xj %in% j)
          nd <- c(length(i), length(j))
          x <- new("dgTMatrix", Dim = nd,
                   i = match(xi[sel], i) - 1:1,
                   j = match(xj[sel], j) - 1:1,
                   x = x@x[sel])
          if (drop && any(nd == 1)) drop(as(x,"matrix")) else x
      })

setMethod("[", signature(x = "dgTMatrix",
			 i = "numeric", j = "numeric", drop = "missing"),
	  function(x,i,j, drop) callGeneric(x,i,j,drop= TRUE))


setMethod("crossprod", signature(x = "dgTMatrix", y = "missing"),
          function(x, y = NULL)
          .Call("csc_crossprod", as(x, "dgCMatrix")))

setMethod("crossprod", signature(x = "dgTMatrix", y = "matrix"),
          function(x, y = NULL)
          .Call("csc_matrix_crossprod", as(x, "dgCMatrix"), y))

##setMethod("crossprod", signature(x = "dgTMatrix", y = "numeric"),
##          function(x, y = NULL)
##          .Call("csc_matrix_crossprod", as(x, "dgCMatrix"), as.matrix(y)))

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
          function(object, ...)
              isTRUE(all.equal(as(object, "dgCMatrix"),
                               as(t(object), "dgCMatrix"))))

setAs("dgTMatrix", "dsCMatrix",
      function(from) {
          if (!isSymmetric(from))
              stop("cannot coerce non-symmetric matrix to dsCMatrix class")
          upper <- from@i <= from@j
          uC <- as(new("dgTMatrix", Dim = from@Dim, i = from@i[upper],
                       j = from@j[upper], x = from@x[upper]), "dgCMatrix")
          new("dsCMatrix", Dim = uC@Dim, p = uC@p, i = uC@i, x = uC@x, uplo = "U")
      })

setAs("matrix", "dgTMatrix",
      function(from) {
          x <- as.double(from)
          nz <- as.logical(x)
          new("dgTMatrix", Dim = dim(from),
              i = as.integer(row(from) - 1)[nz] ,
              j = as.integer(col(from) - 1)[nz],
              x = x[nz])
      })

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
